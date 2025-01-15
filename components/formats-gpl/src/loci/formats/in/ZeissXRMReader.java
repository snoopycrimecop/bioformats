/*
 * #%L
 * OME Bio-Formats package for reading and converting biological file formats.
 * %%
 * Copyright (C) 2005 - 2025 Open Microscopy Environment:
 *   - Board of Regents of the University of Wisconsin-Madison
 *   - Glencoe Software, Inc.
 *   - University of Dundee
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 2 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
 * #L%
 */

package loci.formats.in;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import loci.common.Constants;
import loci.common.DataTools;
import loci.common.DateTools;
import loci.common.Location;
import loci.common.RandomAccessInputStream;
import loci.common.services.DependencyException;
import loci.common.services.ServiceFactory;
import loci.formats.CoreMetadata;
import loci.formats.FormatException;
import loci.formats.FormatReader;
import loci.formats.FormatTools;
import loci.formats.MetadataTools;
import loci.formats.meta.MetadataStore;
import loci.formats.services.POIService;
import loci.formats.services.POIServiceImpl;

/**
 * ZeissXRMReader is the file format reader for Zeiss X-Ray Microscopy
 * .txm and .txrm files.
 */
public class ZeissXRMReader extends FormatReader {

  // -- Constants --

  private static final String IMAGE_INFO = "Root Entry/ImageInfo/";

  // -- Fields --

  private transient POIService poi;
  private List<String> imagePaths = new ArrayList<String>();

  // -- Constructor --

  /** Constructs a new Zeiss XRM reader. */
  public ZeissXRMReader() {
    super("Zeiss XRM", new String[] {"txm", "txrm"});
    domains = new String[] {FormatTools.UNKNOWN_DOMAIN};
    suffixSufficient = true;
  }

  // -- IFormatReader API methods --

  /* @see loci.formats.IFormatReader#isThisType(RandomAccessInputStream) */
  @Override
  public boolean isThisType(RandomAccessInputStream stream) throws IOException {
    final int blockLen = 4;
    if (!FormatTools.validStream(stream, blockLen, false)) return false;
    return stream.readInt() == POIServiceImpl.POI_MAGIC_BYTES;
  }

  /**
   * @see loci.formats.IFormatReader#openBytes(int, byte[], int, int, int, int)
   */
  @Override
  public byte[] openBytes(int no, byte[] buf, int x, int y, int w, int h)
    throws FormatException, IOException
  {
    FormatTools.checkPlaneParameters(this, no, buf.length, x, y, w, h);

    if (poi == null) {
      initPOIService();
    }
    int bpp = FormatTools.getBytesPerPixel(getPixelType());

    try (RandomAccessInputStream stream = poi.getDocumentStream(imagePaths.get(no))) {
      int skipBeginRow = x * bpp;
      int skipEndRow = bpp * (getSizeX() - w - x);
      int rowLen = w * bpp;
      for (int row=h-1; row>=0; row--) {
        stream.skipBytes(skipBeginRow);
        stream.read(buf, row * rowLen, rowLen);
        stream.skipBytes(skipEndRow);
      }
    }

    return buf;
  }

  /* @see loci.formats.IFormatReader#close(boolean) */
  @Override
  public void close(boolean fileOnly) throws IOException {
    super.close(fileOnly);
    if (!fileOnly) {
      if (poi != null) poi.close();
      poi = null;
      imagePaths.clear();
    }
  }

  // -- Internal FormatReader API methods --

  /* @see loci.formats.FormatReader#initFile(String) */
  @Override
  protected void initFile(String id) throws FormatException, IOException {
    super.initFile(id);

    CoreMetadata m = core.get(0);

    initPOIService();

    final List<String> allFiles = poi.getDocumentList();
    if (allFiles.isEmpty()) {
      throw new FormatException(
        "No files were found - the .cxd may be corrupt.");
    }

    for (String name : allFiles) {
      if (name.startsWith("Root Entry/ImageData")) {
        int index = Integer.parseInt(name.substring(name.lastIndexOf("Image") + 5)) - 1;
        while (index >= imagePaths.size()) {
          imagePaths.add(null);
        }
        imagePaths.set(index, name);

        // this is just a raw image with no metadata, so don't try to open it here
        continue;
      }
      try (RandomAccessInputStream stream = poi.getDocumentStream(name)) {
        stream.order(true);

        if (name.equals(IMAGE_INFO + "ImageWidth")) {
          m.sizeX = stream.readInt();
        }
        else if (name.equals(IMAGE_INFO + "ImageHeight")) {
          m.sizeY = stream.readInt();
        }
        else if (name.equals(IMAGE_INFO + "DataType")) {
          m.pixelType = getPixelType(stream.readInt());
        }
        ///* debug */ System.out.println(name + " (" + stream.length() + " bytes)");
      }
    }

    m.sizeZ = imagePaths.size();
    m.sizeT = 1;
    m.sizeC = 1;
    m.imageCount = getSizeZ() * getSizeC() * getSizeT();
    m.dimensionOrder = "XYZTC";
    m.littleEndian = true;

    MetadataStore store = makeFilterMetadata();
    MetadataTools.populatePixels(store, this, true);
  }

  // -- Helper methods --

  private void initPOIService() throws FormatException, IOException {
   try {
      ServiceFactory factory = new ServiceFactory();
      poi = factory.getInstance(POIService.class);
    }
    catch (DependencyException de) {
      throw new FormatException("POI library not found", de);
    }

    poi.initialize(Location.getMappedId(getCurrentFile()));
  }

  private int getPixelType(int dataType) throws FormatException {
    switch (dataType) {
      case 2:
        return FormatTools.INT8;
      case 3:
        return FormatTools.UINT8;
      case 4:
        return FormatTools.INT16;
      case 5:
        return FormatTools.UINT16;
      case 6:
        return FormatTools.INT32;
      case 7:
        return FormatTools.UINT32;
      case 10:
        return FormatTools.FLOAT;
      case 11:
        return FormatTools.DOUBLE;
    }
    throw new FormatException("Unsupported data type: " + dataType);
  }

}
