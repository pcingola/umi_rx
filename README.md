
# `umi_rx`

Add UMI tag to `RX` tag in a BAM file, also add 'MQ' tag

# Install

### Prerequisites

- Java 12

### Compile

Run `script/make.sh`
The JAR file is created at `target/umi_rx-0.1-jar-with-dependencies.jar`

### Running

```
java -Xmx2G UmiRx.jar in.bam out.bam
```

You can use '-' for STDIN / STDOUT (in.bam or out.bam respectively)
