/*
 * #%L
 * OME Bio-Formats package for reading and converting biological file formats.
 * %%
 * Copyright (C) 2024 Open Microscopy Environment:
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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import loci.common.DataTools;
import loci.common.Location;
import loci.formats.CoreMetadata;
import loci.formats.FormatException;
import loci.formats.FormatReader;
import loci.formats.FormatTools;
import loci.formats.MetadataTools;
import loci.formats.WellContainer;
import loci.formats.meta.MetadataStore;

import org.json.JSONException;
import org.json.JSONArray;
import org.json.JSONObject;

/**
 * JDCEReader is the file format reader for Molecular Devices JDCE plates.
 */
public class JDCEReader extends FormatReader {

  // -- Constants --

  // -- Fields --

  private List<JDCEWell> wells = new ArrayList<JDCEWell>();
  private String imageFileCSV = null;
  private transient MinimalTiffReader helper = new MinimalTiffReader();

  // -- Constructor --

  /** Constructs a new JDCE reader. */
  public JDCEReader() {
    super("Molecular Devices JDCE", new String[] {"jdce"});
    suffixSufficient = true;
    domains = new String[] {FormatTools.HCS_DOMAIN};
    hasCompanionFiles = true;
    datasetDescription = "One .jdce (JSON) file with at least one .tif/.tiff file";
  }

  // -- IFormatReader API methods --

  /* @see loci.formats.IFormatReader#isSingleFile(String) */
  @Override
  public boolean isSingleFile(String id) throws FormatException, IOException {
    return false;
  }

  /* @see loci.formats.IFormatReader#fileGroupOption(String) */
  @Override
  public int fileGroupOption(String id) throws FormatException, IOException {
    return FormatTools.MUST_GROUP;
  }

  /**
   * @see loci.formats.IFormatReader#openBytes(int, byte[], int, int, int, int)
   */
  @Override
  public byte[] openBytes(int no, byte[] buf, int x, int y, int w, int h)
    throws FormatException, IOException
  {
    FormatTools.checkPlaneParameters(this, no, buf.length, x, y, w, h);
    Arrays.fill(buf, getFillColor());

    String file = getFile(no);
    helper.setId(file);
    return helper.openBytes(0, buf, x, y, w, h);
  }

  /* @see loci.formats.IFormatReader#getSeriesUsedFiles(boolean) */
  @Override
  public String[] getSeriesUsedFiles(boolean noPixels) {
    FormatTools.assertId(currentId, true, 1);
    // TODO
    return new String[] {currentId, imageFileCSV};
  }

  /* @see loci.formats.IFormatReader#close(boolean) */
  @Override
  public void close(boolean fileOnly) throws IOException {
    super.close(fileOnly);
    if (helper != null) {
      helper.close(fileOnly);
    }
    if (!fileOnly) {
      imageFileCSV = null;
      if (wells != null) {
        wells.clear();
      }
    }
  }

  /* @see loci.formats.IFormatReader#getOptimalTileWidth() */
  @Override
  public int getOptimalTileWidth() {
    FormatTools.assertId(currentId, true, 1);
    // TODO
    return super.getOptimalTileWidth();
  }

  /* @see loci.formats.IFormatReader#getOptimalTileHeight() */
  @Override
  public int getOptimalTileHeight() {
    FormatTools.assertId(currentId, true, 1);
    // TODO
    return super.getOptimalTileHeight();
  }

  // -- Internal FormatReader API methods --

  /* @see loci.formats.FormatReader#initFile(String) */
  @Override
  protected void initFile(String id) throws FormatException, IOException {
    super.initFile(id);

    CoreMetadata ms0 = core.get(0);
    try {
      JSONObject root = new JSONObject(DataTools.readFile(id));

      JSONObject imageStack = root.getJSONObject("ImageStack");
      if (imageStack == null) {
        throw new FormatException("Could not find image stack definition");
      }

      String imageFormat = imageStack.getString("ImageFormat");
      if (!"TIFF".equalsIgnoreCase(imageFormat)) {
        throw new FormatException("Unsupported image format " + imageFormat);
      }

      JSONObject acquisition = imageStack.getJSONObject("AutoLeadAcquisitionProtocol");
      if (acquisition == null) {
        throw new FormatException("Could not find acquisition definition");
      }

      JSONObject objective = acquisition.getJSONObject("ObjectiveCalibration");
      // TODO: parse objective

      JSONObject plate = acquisition.getJSONObject("Plate");
      // TODO: parse plate

      JSONObject plateMap = acquisition.getJSONObject("PlateMap");
      if (plateMap == null) {
        throw new FormatException("Could not find plate map, cannot determine dimensions");
      }
      JSONObject timeSchedule = plateMap.getJSONObject("TimeSchedule");
      if (timeSchedule == null) {
        throw new FormatException("Could not find time schedule, cannot determine SizeT");
      }
      ms0.sizeT = timeSchedule.getInt("NumberOfTimepoints");

      JSONObject zDimension = plateMap.getJSONObject("ZDimensionParameters");
      if (zDimension == null) {
        throw new FormatException("Could not find Z dimension parameters, cannot determine SizeZ");
      }
      ms0.sizeZ = zDimension.getInt("NumberOfSlices");

      JSONArray wavelengths = acquisition.getJSONArray("Wavelengths");
      if (wavelengths == null) {
        throw new FormatException("Could not find wavelength array, cannot determine SizeC");
      }
      ms0.sizeC = wavelengths.length();

      JSONArray metadataFiles = imageStack.getJSONArray("ImageMetadataFiles");
      if (metadataFiles == null || metadataFiles.length() == 0) {
        throw new FormatException("Could not find image metadata CSV, cannot get list of TIFF files");
      }
      imageFileCSV = metadataFiles.getString(0);
    }
    catch (JSONException e) {
      throw new FormatException("Could not parse .jdce file", e);
    }
    ms0.imageCount = getSizeZ() * getSizeC() * getSizeT();
    ms0.dimensionOrder = "XYCZT";

    if (imageFileCSV == null) {
      throw new FormatException("Image metadata CSV not found, cannot get list of TIFF files");
    }

    String[] csvLines = DataTools.readFile(imageFileCSV).split("\r\n");
    List<String> columns = Arrays.asList(csvLines[0].split(","));
    int wellRowIndex = columns.indexOf("Row");
    int wellColIndex = columns.indexOf("Column");
    int fieldIndex = columns.indexOf("Field");
    int wavelengthIndex = columns.indexOf("Wavelength");
    int timepointIndex = columns.indexOf("Timepoint");
    int zIndex = columns.indexOf("ZIndex");
    int subfolderIndex = columns.indexOf("ImageSubFolderPath");
    int fileNameIndex = columns.indexOf("ImageFileName");

    Location parentDir = new Location(getCurrentFile()).getAbsoluteFile().getParentFile();
    JDCEWell currentWell = null;
    boolean firstFile = true;
    for (int i=1; i<csvLines.length; i++) {
      String[] line = csvLines[i].split(",");

      int[] position = new int[4];
      position[0] = Integer.parseInt(line[fieldIndex]);
      position[1] = Integer.parseInt(line[zIndex]);
      position[2] = Integer.parseInt(line[wavelengthIndex]);
      position[3] = Integer.parseInt(line[timepointIndex]);

      String subfolder = line[subfolderIndex];
      String filename = line[fileNameIndex];
      Location subfolderFile = new Location(parentDir, subfolder);
      String imagePath = new Location(subfolderFile, filename).getAbsolutePath();

      int wellRow = Integer.parseInt(line[wellRowIndex]) - 1;
      int wellCol = Integer.parseInt(line[wellColIndex]) - 1;
      if (currentWell == null || currentWell.getRowIndex() != wellRow ||
        currentWell.getColumnIndex() != wellCol)
      {
        currentWell = lookupWell(wellRow, wellCol);
      }
      currentWell.addFile(imagePath, position);
      currentWell.setFieldCount((int) Math.max(currentWell.getFieldCount(), position[0] + 1));

      if (firstFile) {
        try {
          helper.setId(imagePath);
          CoreMetadata m = helper.getCoreMetadataList().get(0);
          ms0.sizeX = m.sizeX;
          ms0.sizeY = m.sizeY;
          ms0.pixelType = m.pixelType;
          ms0.littleEndian = m.littleEndian;
          ms0.sizeC *= m.sizeC;
          ms0.rgb = m.rgb;

          firstFile = false;
        }
        catch (FormatException | IOException e) {
          LOGGER.debug("Could not read " + imagePath, e);
        }
      }
    }
    for (JDCEWell well : wells) {
      for (int f=0; f<well.getFieldCount(); f++) {
        core.add(new CoreMetadata(ms0));
      }
    }
    core.remove(0);

    MetadataStore store = makeFilterMetadata();
    MetadataTools.populatePixels(store, this, true);

    wells.sort(null);

    int imageIndex = 0;
    for (int w=0; w<wells.size(); w++) {
      wells.get(w).fillMetadataStore(store, 0, 0, w, 0, imageIndex);
      imageIndex += wells.get(w).getFieldCount();
    }
  }

  private JDCEWell lookupWell(int row, int col) {
    for (JDCEWell well : wells) {
      if (well.getRowIndex() == row && well.getColumnIndex() == col) {
        return well;
      }
    }
    JDCEWell well = new JDCEWell(row, col);
    wells.add(well);
    return well;
  }

  private String getFile(int no) {
    int index = 0;
    for (JDCEWell well : wells) {
      if (well.getFieldCount() + index > getSeries()) {
        return well.getFile(getSeries() - index, no);
      }
      index += well.getFieldCount();
    }
    return null;
  }

  class JDCEWell extends WellContainer {
    HashMap<int[], Integer> posMap = new HashMap<int[], Integer>();

    public JDCEWell(int row, int col) {
      super(0, row, col, 1);
    }

    /**
     * Add file at position [field, z, c, t].
     */
    public void addFile(String file, int[] position) {
      super.addFile(file);
      posMap.put(position, getAllFiles().size() - 1);
    }

    public String[] getFiles(int fieldIndex) {
      List<String> allFiles = getAllFiles();
      String[] fieldFiles = new String[getImageCount()];
      for (int[] key : posMap.keySet()) {
        if (key[0] == fieldIndex) {
          fieldFiles[getIndex(key[1], key[2], key[3])] = allFiles.get(posMap.get(key));
        }
      }
      return fieldFiles;
    }
  }
}
