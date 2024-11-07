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
import java.util.List;

import loci.common.DataTools;
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

  private List<WellContainer> wells = new ArrayList<WellContainer>();
  private String imageFileCSV = null;

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

    return buf;
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

    if (imageFileCSV == null) {
      throw new FormatException("Image metadata CSV not found, cannot get list of TIFF files");
    }

    String[] csvLines = DataTools.readFile(imageFileCSV).split("\r\n");
    String[] columns = csvLines[0].split(",");
    for (int i=1; i<csvLines.length; i++) {
      String[] line = csvLines[i].split(",");
    }

    MetadataStore store = makeFilterMetadata();
    MetadataTools.populatePixels(store, this, true);
  }
}
