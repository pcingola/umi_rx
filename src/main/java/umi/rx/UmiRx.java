package umi.rx;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.BlockCompressedOutputStream;

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
	public static long MAX_READS = 1000000;
	public static long SHOW_EVERY = 1000000;

	boolean calcRx, calcMc, calcMq;
	int compressionLevel;
	long countMc, countMq, countRx;
	boolean debug = false;
	String inBam, outBam;
	SamReader samReader;
	SAMFileWriter samWriter;
	boolean useSamOutput;
	boolean verbose = true;

	public static void main(String[] args) {
		// Parse command line
		CommandLine cmd = parseCommandLine(args);
		String in = cmd.hasOption('i') ? cmd.getOptionValue('i') : "-";
		String out = cmd.hasOption('o') ? cmd.getOptionValue('o') : "-";

		// Process input BAM, write output BAM
		UmiRx umirx = new UmiRx(in, out);

		// Compression level
		if (cmd.hasOption("comp")) {
			String compStr = cmd.getOptionValue("comp");
			int level = Integer.parseInt(compStr);
			umirx.setCompressionLevel(level);
		}

		// Other options
		umirx.setDebug(cmd.hasOption('v'));
		umirx.setVerbose(cmd.hasOption('d'));
		umirx.setUseSamOutput(cmd.hasOption('s'));
		umirx.setCalcMc(cmd.hasOption('c'));
		umirx.setCalcMq(cmd.hasOption('q'));
		umirx.setCalcRx(cmd.hasOption('x'));

		// Process BAM
		umirx.open();
		umirx.transform();
		umirx.close();
	}

	/**
	 * Parse command line options
	 */
	protected static CommandLine parseCommandLine(String[] args) {
		// Create command line options
		Options options = new Options();
		options.addOption(new Option("h", "help", false, "Help"));
		options.addOption(new Option("d", "debug", false, "Debug mode, only porcess " + MAX_READS + " reads"));
		options.addOption(new Option("v", "verbose", false, "Verbose"));
		options.addOption(new Option("c", "mc", false, "Add MC tag (mate cigar)"));
		options.addOption(new Option("q", "mq", false, "Add MQ tag (mate mapping quality)"));
		options.addOption(new Option("x", "rx", false, "Add RX tag (UMI from read name)"));
		options.addOption(new Option("s", "sam", false, "Use SAM output format instead ob BAM"));

		Option comp = new Option("l", "comp", true, "Compression level for output BAM");
		comp.setArgName("level");
		options.addOption(comp);

		Option in = new Option("i", "in", true, "Input BAM. Default: STDIN");
		comp.setArgName("level");
		options.addOption(in);

		Option out = new Option("o", "out", true, "Output BAM/SAM. Default: STDIN");
		comp.setArgName("level");
		options.addOption(out);

		// Parse command line options
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = null;
		boolean showHelp = false;

		// Parse args, show help if needed or if there are errors parsing
		try {
			cmd = parser.parse(options, args);

			// Check command line options
			if (cmd.hasOption('h')) {
				// Show command help and exit
				showHelp = true;
			} else if (cmd.getArgList().size() != 0) {
				System.out.println("Error: Unknown parameters " + cmd.getArgList());
				showHelp = true;
			} else if (!(cmd.hasOption('c') || cmd.hasOption('q') || cmd.hasOption('x'))) {
				System.out.println("Error: At least one of the options {'--mc', '--mq', '--rx'} must be specified");
				showHelp = true;
			}
		} catch (ParseException e) {
			showHelp = true;
			System.out.println(e.getMessage());
		} finally {
			if (showHelp) {
				HelpFormatter formatter = new HelpFormatter();
				formatter.printHelp("java -jar UmiRx.jar [OPTIONS]", options);
				System.exit(1);
			}
		}

		return cmd;
	}

	public UmiRx(String inBam, String outBam) {
		this.inBam = inBam;
		this.outBam = outBam;
		compressionLevel = BlockCompressedOutputStream.getDefaultCompressionLevel();
		countMc = countMq = countRx = 0;
	}

	/**
	 * Add MC tag, if it doesn't exits
	 */
	void addMc(SAMRecord sr, String cigarMate) {
		if (calcMc && sr.getAttribute(MC) == null) {
			sr.setAttribute(MC, cigarMate);
			countMc++;
		}
	}

	/**
	 * Add MQ tag, if it doesn't exits
	 */
	void addMq(SAMRecord sr, int qualMate) {
		if (calcMq && sr.getAttribute(MQ) == null) {
			sr.setAttribute(MQ, qualMate);
			countMq++;
		}
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
		if (!calcRx || sr.getAttribute(RX) != null) return;
		String readName = sr.getReadName();
		int idx = readName.lastIndexOf(':');
		if (idx < 0) throw new RuntimeException("Could not find ':' in read name. Read: " + sr);
		String umi = readName.substring(idx + 1);
		sr.setAttribute(RX, umi);
		countRx++;
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

	public int getCompressionLevel() {
		return compressionLevel;
	}

	public boolean isDebug() {
		return debug;
	}

	public boolean isUseSamOutput() {
		return useSamOutput;
	}

	public boolean isVerbose() {
		return verbose;
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
		SAMFileWriterFactory swf = new SAMFileWriterFactory();
		if (useSamOutput) {
			samWriter = outBam.equals("-") ? swf.makeSAMWriter(samHeader, false, System.out) : swf.makeSAMWriter(samHeader, false, new File(outBam));
		} else {
			samWriter = outBam.equals("-") ? swf.makeBAMWriter(samHeader, false, System.out) : swf.makeBAMWriter(samHeader, false, new File(outBam), compressionLevel);
		}
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
		addMc(sr, "");
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
			} else {
				System.err.println("WARNIGN: Neither first nor second pair " + sr);
				addMq(sr, 0);
				addMc(sr, "");
			}

			// Save
			samWriter.addAlignment(sr);
		}
	}

	public void setCalcMc(boolean calcMc) {
		this.calcMc = calcMc;
	}

	public void setCalcMq(boolean calcMq) {
		this.calcMq = calcMq;
	}

	public void setCalcRx(boolean calcRx) {
		this.calcRx = calcRx;
	}

	public void setCompressionLevel(int compressionLevel) {
		this.compressionLevel = compressionLevel;
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public void setUseSamOutput(boolean useSamOutput) {
		this.useSamOutput = useSamOutput;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	public void transform() {
		if (calcMc || calcMq) transformM();
		else transformRx();
	}

	/**
	 * Transform records from inBam and write them to outBam
	 * Add MC/MQ tags.
	 * Note: We need to analyze all reads with the same read name at the same time
	 */
	protected void transformM() {
		long readNum = 0;
		List<SAMRecord> srs = new ArrayList<>();

		String readNamePrev = "";
		for (SAMRecord sr : samReader) {

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
			if (readNum % SHOW_EVERY == 0) System.err.println(readNum + "\t" + sr + "\tcountMc: " + countMc + "\tcountMq: " + countMq + "\tcountRx: " + countRx);
			if (debug && readNum > MAX_READS) {
				System.err.println("WARNING: Debug mode, breaking after " + MAX_READS + " reads");
				break;
			}
			readNum++;
		}

		process(srs); // Process last list of reads
		System.err.println(readNum + "\tcountMc: " + countMc + "\tcountMq: " + countMq + "\tcountRx: " + countRx);
	}

	/**
	 * Only add RX tag
	 * Note: We only need to analyze one read at a time
	 */
	protected void transformRx() {
		long readNum = 0;
		for (SAMRecord sr : samReader) {

			// Add RX tag and write read
			addRx(sr);
			samWriter.addAlignment(sr);

			// Show progress
			if (readNum % SHOW_EVERY == 0) System.err.println(readNum + "\t" + sr + "\tcountMc: " + countMc + "\tcountMq: " + countMq + "\tcountRx: " + countRx);
			if (debug && readNum > MAX_READS) {
				System.err.println("WARNING: Debug mode, breaking after " + MAX_READS + " reads");
				break;
			}
			readNum++;
		}
		System.err.println(readNum + "\tcountMc: " + countMc + "\tcountMq: " + countMq + "\tcountRx: " + countRx);
	}

}
