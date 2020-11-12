package umi.rx;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

/**
 * Some custom processing needed for some BAM files
 * particularly when trying to process Dragen generated BAMs in FGBIO.
 *
 * IMPORTANT: Reads in the input BAM MUST be sorted by read-name (not chr:pos)
 *
 *
 * - Add UMIs from read name into an 'RX' tag
 * - Add MQ tag
 *
 * @author pcingola
 *
 */
public class UmiRx {

	public static final String MC = SAMTag.MC.name();
	public static final String MQ = SAMTag.MQ.name();
	public static final String RX = SAMTag.RX.name();

	public static final String HOME = System.getProperty("user.home");
	public static long MAX_READS = 100000;
	public static long SHOW_EVERY = 100000;

	boolean debug = false;
	boolean verbose = true;
	String inBam, outBam;
	SamReader samReader;
	SAMFileWriter samWriter;

	public static void main(String[] args) {
		// Parse command line options
		if (args.length != 2) {
			System.err.println("Usage: java -jar UmiRx.jar input.bam output.bam");
			System.exit(1);
		}
		String in = args[0];
		String out = args[1];

		// Process input BAM, write output BAM
		UmiRx umirx = new UmiRx(in, out);
		umirx.open();
		umirx.transform();
		umirx.close();
	}

	public UmiRx(String inBam, String outBam) {
		this.inBam = inBam;
		this.outBam = outBam;
	}

	/**
	 * Add MC tag, if it doesn't exits
	 */
	void addMc(SAMRecord sr, String cigarMate) {
		if (sr.getAttribute(MC) == null) sr.setAttribute(MC, cigarMate);
	}

	/**
	 * Add MQ tag, if it doesn't exits
	 */
	void addMq(SAMRecord sr, int qualMate) {
		if (sr.getAttribute(MQ) == null) sr.setAttribute(MQ, qualMate);
	}

	/**
	 * Add RX tag, if it doesn't exits.
	 *
	 * Get UMI from read name and add it as RX tag
	 *
	 * UMI is the last entry in the read name (when splitting by ':')
	 * Example:
	 * 		Read name: A00324:79:HJ5CMDSXX:2:1101:19705:1172:CGCACG
	 *      UMI      : CGCACG
	 */
	protected void addRx(SAMRecord sr) {
		if (sr.getAttribute(RX) != null) return;
		String readName = sr.getReadName();
		int idx = readName.lastIndexOf(':');
		if (idx < 0) throw new RuntimeException("Could not find ':' in read name. Read: " + sr);
		String umi = readName.substring(idx + 1);
		sr.setAttribute(RX, umi);
	}

	/**
	 * Close SAM reader and writer
	 */
	public void close() {
		try {
			samWriter.close();
			samReader.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Open input and output BAMs
	 */
	public void open() {
		// Open input BAM
		SamInputResource samIn = inBam.equals("-") ? SamInputResource.of(System.in) : SamInputResource.of(inBam);

		samReader = SamReaderFactory.make().validationStringency(ValidationStringency.LENIENT).open(samIn);
		SAMFileHeader samHeader = samReader.getFileHeader();

		// Create output BAM file
		//		SamInputResource samIn = inBam.equals("-") ? SamInputResource.of(System.in) : SamInputResource.of(inBam);
		SAMFileWriterFactory swf = new SAMFileWriterFactory();
		samWriter = outBam.equals("-") ? swf.makeBAMWriter(samHeader, false, System.out) : swf.makeBAMWriter(samHeader, false, new File(outBam));
	}

	protected void process(List<SAMRecord> srs) {
		switch (srs.size()) {
		case 0:
			// No reads, nothing to do
			break;
		case 1:
			process1(srs);
			break;
		case 2:
			process2(srs);
			break;
		default:
			process3orMore(srs);
		}
	}

	/**
	 * Process list of reads containing only one read
	 */
	protected void process1(List<SAMRecord> srs) {
		SAMRecord sr = srs.get(0);
		addRx(sr);
		addMq(sr, 0);
		samWriter.addAlignment(sr);
	}

	/**
	 * Process list of reads containing two reads
	 * Note: This is the most common case (assuming pair-end reads)
	 */
	protected void process2(List<SAMRecord> srs) {
		SAMRecord sr1 = srs.get(0);
		SAMRecord sr2 = srs.get(1);

		int mq1 = sr1.getMappingQuality();
		int mq2 = sr2.getMappingQuality();
		String cigar1 = sr1.getCigarString();
		String cigar2 = sr2.getCigarString();

		// MQ tags is the mapping quality of "mate read"
		addMq(sr1, mq2);
		addMq(sr2, mq1);

		// MC tags is the cigar of "mate read"
		addMc(sr1, cigar2);
		addMc(sr2, cigar1);

		// Add UMIs
		addRx(sr1);
		addRx(sr2);

		// Save
		samWriter.addAlignment(sr1);
		samWriter.addAlignment(sr2);
	}

	/**
	 * Process a pair of SAM records
	 */
	protected void process3orMore(List<SAMRecord> srs) {
		boolean ok = true;
		int mq1 = 0, mq2 = 0;
		String cigar1 = "", cigar2 = "";

		// Get mapping qualities
		for (SAMRecord sr : srs) {
			if (sr.getReadUnmappedFlag() || sr.getMateUnmappedFlag()) {
				ok = false;
			} else if (sr.getFirstOfPairFlag()) {
				mq1 = Math.max(mq1, sr.getMappingQuality());
				if (cigar1.isEmpty() || !sr.isSecondaryOrSupplementary()) cigar1 = sr.getCigarString();
			} else if (sr.getSecondOfPairFlag()) {
				mq2 = Math.max(mq1, sr.getMappingQuality());
				if (cigar2.isEmpty() || !sr.isSecondaryOrSupplementary()) cigar2 = sr.getCigarString();
			}
		}

		if (!ok) mq1 = mq2 = 0;

		// Set RX, MQ and save
		for (SAMRecord sr : srs) {
			addRx(sr); // Add RX tag

			// Add MQ tag (mapping quality of paired read)
			if (sr.getFirstOfPairFlag()) {
				addMq(sr, mq2);
				addMc(sr, cigar2);
			} else if (sr.getSecondOfPairFlag()) {
				addMq(sr, mq1);
				addMc(sr, cigar1);
			}

			// Save
			samWriter.addAlignment(sr);
		}
	}

	/**
	 * Transform all records from inBam and write them to outBam
	 */
	public void transform() {
		long readNum = 1;
		List<SAMRecord> srs = new ArrayList<>();

		String readNamePrev = "";
		for (SAMRecord sr : samReader) {
			readNum++;

			// Collect all reads with the same name in a list
			// Process the list of reads when the read name changes
			String readName = sr.getReadName();
			if (!readName.equals(readNamePrev)) {
				// Read name changed, process reads, then clear list
				process(srs);
				srs.clear();
			}
			srs.add(sr);

			// Prepare for next iteration
			readNamePrev = readName;

			// Show progress
			if (readNum % SHOW_EVERY == 0) System.err.println(readNum + "\t" + sr);
			if (debug && readNum > MAX_READS) {
				System.err.println("WARNING: Debug mode, breaking after " + MAX_READS + " reads");
				break;
			}
		}

		process(srs); // Process last list of reads

	}

}
