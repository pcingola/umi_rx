#include <stdio.h>
#include <unistd.h>

#include "htslib/sam.h"
#include "htslib/vcf.h"

int main(int argc, char **argv) {
    if(argc != 3) {
        fprintf(stderr, "Usage: %s input.bam output.bam\n", argv[0]);
        return 1;
    }

    char moder[8];
    char *filein = argv[1];
    char *fileout = argv[2];

    // Open in.bam
    htsFile *in = hts_open(filein, "r");
    if (!in) {
        fprintf(stderr, "Error opening \"%s\"\n", filein);
        exit(1);
    }

    // Open out.bam
    htsFile *out = hts_open(fileout, "wb");
    if (!out) {
        fprintf(stderr, "Error opening \"%s\"\n", filein);
        exit(1);
    }

    // Read header
    sam_hdr_t *header = sam_hdr_read(in);
    if (header == NULL) {
        fprintf(stderr, "Couldn't read header for \"%s\"\n", filein);
        exit(1);
    }

    // Write header
    if (sam_hdr_write(out, header) < 0) {
        fprintf(stderr, "Error writing output header.\n");
        exit(1);
    }

    bam1_t *aln = bam_init1();
    long read_num;
    for (read_num = 0; sam_read1(in, header, aln) > 0; read_num++) {
        char *read_name = bam_get_qname(aln);
        int32_t pos = aln->core.pos + 1;
        char *chr = header->target_name[aln->core.tid];

        // Find UMI part
        char *umi = strrchr(read_name, ':');
        if (!umi) {
            fprintf(stderr, "Error: Could not find UMI from read name, read_number=%ld, chr='%s', pos=%d, read_name='%s'\n", read_num, chr, pos, read_name);
            exit(1);
        }
        umi++; // We want the string starting right after ':'

        // Check UMI length
        int umilen = strlen(umi) + 1;
        if(umilen != 7) {
            fprintf(stderr, "WARNING: UMI len is %d, read name, read_number=%ld, chr='%s', pos=%d, read_name='%s', umi(length=%d)='%s' \n", umilen, read_num, chr, pos, read_name, umilen, umi);
        }

        // Show every N reads
        if( read_num % 1 == 0 ) {
            printf("read_number=%ld, chr='%s', pos=%d, read_name='%s', umi(length=%d)='%s'\n", read_num, chr, pos, read_name, umilen, umi);
        }

        // Add UMI to 'RX' tag
        printf("bam_aux_append\tread_number=%ld'\n", read_num);
        if (bam_aux_append(aln, "RX", 'Z', umilen, (uint8_t *) umi) < 0) {
            fprintf(stderr, "Error updating RX tag");
            exit(1);
        }

        // Write alignment to output
        if (sam_write1(out, header, aln) < 0) {
            fprintf(stderr, "Error writing output alignment, read_number=%ld, chr='%s', pos=%d, read_name='%s'\n", read_num, chr, pos, read_name);
            exit(1);
        }
    }

    printf("Finished: %ld reads processed\n", read_num);

    // Close files
    if (hts_close(out) < 0) {
        fprintf(stderr, "Error closing \"%s\"\n", fileout);
        exit(1);
    }

    if (hts_close(in) < 0) {
        fprintf(stderr, "Error closing \"%s\"\n", filein);
        exit(1);
    }

    // Free memory
    bam_destroy1(aln);
    sam_hdr_destroy(header);

    return 0;
}
