package loci.formats.filter;

import ucar.nc2.filter.FilterProvider;
import ucar.nc2.filter.Filter;

import io.airlift.compress.lz4.Lz4Decompressor;

import java.util.Map;
import java.io.IOException;
import java.nio.ByteBuffer;

public class LZ4Filter extends Filter {

  static final String name = "LZ4 filter";
  static final int id = 32004;
  private Lz4Decompressor decompressor;

  public LZ4Filter(Map<String, Object> properties) {
    decompressor = new Lz4Decompressor();
  }

  @Override
  public String getName() {
    return name;
  }

  @Override
  public int getId() {
    return id;
  }

  @Override
  public byte[] encode(byte[] dataIn) throws IOException {
    // only decoding implementation needed for reading
    byte[] dataOut = dataIn;
    return dataOut;
  }

  @Override
  public byte[] decode(byte[] dataIn) throws IOException {
    ByteBuffer byteBuffer = ByteBuffer.wrap(dataIn);
    
    long totalDecompressedSize = byteBuffer.getLong();
    int decompressedBlockSize = byteBuffer.getInt();
    int compressedBlockSize = byteBuffer.getInt();

    byte[] decompressedBlock = new byte[Math.toIntExact(totalDecompressedSize)];
    byte[] compressedBlock = new byte[compressedBlockSize];

    byteBuffer.get(compressedBlock, 0, compressedBlockSize);
    decompressor.decompress(compressedBlock, 0, compressedBlockSize, decompressedBlock, 0, decompressedBlockSize);

    return decompressedBlock;
  }

  public static class LZ4FilterProvider implements FilterProvider {

    @Override
    public String getName() {
      return name;
    }

    @Override
    public int getId() {
      return id;
    }

    @Override
    public Filter create(Map<String, Object> properties) {
      return new LZ4Filter(properties);
    }

  }

}
