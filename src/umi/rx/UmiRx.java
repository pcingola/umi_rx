package umi.rx;

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

/**
 * Some custom processing needed for some BAM files
 * particularly when trying to process Dragen generated BAMs in FGBIO.
 *
 * - Add UMIs from read name into an 'RX' tag
 * - Add MQ tag
 *
 * @author pcingola
 *
 */
public class UmiRx {

	public static final String HOME = System.getProperty("user.home");
	public static long MAX_READS = 1000;
	public static long SHOW_EVERY = 1000;

	boolean debug = true;
	String inBam, outBam;
	SamReader samReader;
	SAMFileWriter samWriter;

	public static void main(String[] args) {
		String in = HOME + "/umi_rx/in.bam";
		String out = HOME + "/umi_rx/out.bam";

		UmiRx umirx = new UmiRx(in, out);
		System.err.println("Start: Reading " + in);
		umirx.open();
		umirx.transform();
		umirx.close();
		System.err.println("Done");
	}

	public UmiRx(String inBam, String outBam) {
		this.inBam = inBam;
		this.outBam = outBam;
	}

	public void close() {
		try {
			samWriter.close();
			samReader.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	public void open() {
		// Open input BAM
		SamInputResource samIn = SamInputResource.of(inBam);
		samReader = SamReaderFactory.make().validationStringency(ValidationStringency.LENIENT).open(samIn);
		SAMFileHeader samHeader = samReader.getFileHeader();

		// Create output BAM file
		samWriter = (new SAMFileWriterFactory()).makeBAMWriter(samHeader, false, new File(outBam));
	}

	/**
	 * Process a SAM record, transform and write to outBam
	 */
	protected void process(SAMRecord sr) {
		String umi = umi(sr);
		sr.setAttribute("RX", umi);
		System.out.println(umi + "\t" + sr);

		samWriter.addAlignment(sr);
	}

	/**
	 * Transform all records from inBam and write them to outBam
	 */
	public void transform() {
		long readNum = 1;
		for (SAMRecord sr : samReader) {
			process(sr);
			readNum++;
			if (readNum % SHOW_EVERY == 0) System.err.println(sr);
			if (debug && readNum > MAX_READS) {
				System.err.println("WARNING: Debug mode, breaking after " + MAX_READS + " reads");
				break;
			}
		}
	}

	/**
	 * UMI is the last entry in the read name (when splitting by ':')
	 * Example:
	 * 		Read name: A00324:79:HJ5CMDSXX:2:1101:19705:1172:CGCACG
	 *      UMI      : CGCACG
	 */
	protected String umi(SAMRecord sr) {
		String readName = sr.getReadName();
		int idx = readName.lastIndexOf(':');
		if (idx < 0) throw new RuntimeException("Could not find ':' in read name. Read: " + sr);
		return readName.substring(idx + 1);
	}

}
