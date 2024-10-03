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

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import loci.common.Constants;
import loci.common.DataTools;
import loci.common.DateTools;
import loci.common.Location;
import loci.common.Region;
import loci.formats.CoreMetadata;
import loci.formats.FormatException;
import loci.formats.FormatReader;
import loci.formats.FormatTools;
import loci.formats.MetadataTools;
import loci.formats.UnsupportedCompressionException;
import loci.formats.codec.Codec;
import loci.formats.codec.CodecOptions;
import loci.formats.codec.JPEGCodec;
import loci.formats.codec.JPEGXRCodec;
import loci.formats.codec.PassthroughCodec;
import loci.formats.meta.MetadataStore;

import ome.units.UNITS;
import ome.units.quantity.Length;
import ome.units.quantity.Time;
import ome.xml.model.primitives.NonNegativeInteger;
import ome.xml.model.primitives.PositiveInteger;
import ome.xml.model.primitives.Timestamp;

import org.json.JSONException;
import org.json.JSONObject;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.sqlite.SQLiteConfig;

/**
 *
 */
public class TissuegnosticsReader extends FormatReader {

  private static final Logger LOGGER =
    LoggerFactory.getLogger(TissuegnosticsReader.class);

  private List<String> pixelsFiles = new ArrayList<String>();
  private int[] fullResolutionCoreIndex;
  private int[] tileSizeX;
  private int[] tileSizeY;
  private int[] overlapX;
  private int[] overlapY;
  private int[] scaleFactor;
  private int[][] tileRangeX;
  private int[][] tileRangeY;
  private Integer[][] zSteps;

  // -- Constructor --

  /** Constructs a new Tissuegnostics reader.*/
  public TissuegnosticsReader() {
    super("Tissuegnostics", new String[] {"aqproj"});
    hasCompanionFiles = false;
    domains = new String[] {FormatTools.HCS_DOMAIN};
    datasetDescription = "An .aqproj file with one or more .tfcyto database files";
    suffixSufficient = true;
  }

  // -- IFormatReader API methods --

  /* @see loci.formats.IFormatReader#close(boolean) */
  @Override
  public void close(boolean fileOnly) throws IOException {
    super.close(fileOnly);
    if (!fileOnly) {
      pixelsFiles.clear();
      fullResolutionCoreIndex = null;
      tileSizeX = null;
      tileSizeY = null;
      overlapX = null;
      overlapY = null;
      scaleFactor = null;
      tileRangeX = null;
      tileRangeY = null;
      zSteps = null;
    }
  }

  /* @see loci.formats.IFormatReader#getOptimalTileWidth() */
  @Override
  public int getOptimalTileWidth() {
    FormatTools.assertId(currentId, true, 1);
    int region = getRegionIndex();
    return tileSizeX[region] / (int) Math.pow(scaleFactor[region], getLevel());
  }

  /* @see loci.formats.IFormatReader#getOptimalTileHeight() */
  @Override
  public int getOptimalTileHeight() {
    FormatTools.assertId(currentId, true, 1);
    int region = getRegionIndex();
    return tileSizeY[region] / (int) Math.pow(scaleFactor[region], getLevel());
  }

  /* @see loci.formats.IFormatReader#getSeriesUsedFiles(boolean) */
  @Override
  public String[] getSeriesUsedFiles(boolean noPixels) {
    FormatTools.assertId(currentId, true, 1);

    ArrayList<String> files = new ArrayList<String>();
    files.add(getCurrentFile());
    if (!noPixels) {
      files.add(getDBFile());
    }
    return files.toArray(new String[files.size()]);
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

    int regionIndex = getRegionIndex();
    int bpp = FormatTools.getBytesPerPixel(getPixelType());
    int level = getLevel();
    int scale = (int) Math.pow(scaleFactor[regionIndex], level);

    Region dest = new Region(x, y, w, h);
    Connection conn = openConnection(getDBFile());
    try {
      PreparedStatement tiles = conn.prepareStatement(
        "SELECT data, compression, row, column FROM images WHERE level=? AND row>=? AND row<=? AND column>=? AND column<=? AND channel=? AND is_zstack=? AND z_position=? ORDER BY row,column DESC"
      );
      // TODO: account for possiblity of overlap

      int rowOffset = (int) Math.floor((double) tileRangeY[0][regionIndex] / scale);
      int colOffset = (int) Math.floor((double) tileRangeX[0][regionIndex] / scale);

      int startRow = (int) Math.floor((double) y / tileSizeY[regionIndex]) + rowOffset;
      int startCol = (int) Math.floor((double) x / tileSizeX[regionIndex]) + colOffset;
      int endRow = (int) Math.ceil((double) (y + h) / tileSizeY[regionIndex]) + rowOffset;
      int endCol = (int) Math.ceil((double) (x + w) / tileSizeX[regionIndex]) + colOffset;

      tiles.setInt(1, level);
      tiles.setInt(2, startRow);
      tiles.setInt(3, endRow);
      tiles.setInt(4, startCol);
      tiles.setInt(5, endCol);

      int[] zct = getZCTCoords(no);
      tiles.setInt(6, zct[1] + 1);
      tiles.setInt(7, zct[0] > 0 ? 1 : 0);
      tiles.setInt(8, zSteps[regionIndex][zct[0]]);

      // upper left corner of full resolution,
      // in pixels relative to full canvas
      int fullResUpperLeftX = tileRangeX[0][regionIndex] * tileSizeX[regionIndex];
      int fullResUpperLeftY = tileRangeY[0][regionIndex] * tileSizeY[regionIndex];

      // upper left corner of current resolution,
      // in pixels relative to current resolution canvas
      int currentResUpperLeftX = fullResUpperLeftX / scale;
      int currentResUpperLeftY = fullResUpperLeftY / scale;

      int relativeUpperLeftX = currentResUpperLeftX % tileSizeX[regionIndex];
      int relativeUpperLeftY = currentResUpperLeftY % tileSizeY[regionIndex];

      ResultSet subsetTiles = tiles.executeQuery();
      while (subsetTiles.next()) {
        byte[] data = subsetTiles.getBytes(1);
        int compression = subsetTiles.getInt(2);
        int row = subsetTiles.getInt(3);
        int column = subsetTiles.getInt(4);

        LOGGER.debug("found {} bytes for row={}, column={}",
          data.length, row, column);

        CodecOptions options = new CodecOptions();
        options.bitsPerSample = bpp * 8;
        options.width = tileSizeX[regionIndex];
        options.height = tileSizeY[regionIndex];

        Codec codec = getCodec(compression);
        byte[] tile = codec.decompress(data, options);

        int relativeColumn = column - (int) Math.floor((double) tileRangeX[0][regionIndex] / scale);
        int relativeRow = row - (int) Math.floor((double) tileRangeY[0][regionIndex] / scale);
        relativeColumn *= options.width;
        relativeRow *= options.height;

        if (level > 0) {
          relativeRow -= relativeUpperLeftY;
          relativeColumn -= relativeUpperLeftX;
        }

        Region src = new Region(relativeColumn, relativeRow, options.width, options.height);
        Region intersection = src.intersection(dest);

        int outputRowLen = w * bpp;
        int intersectionX = (int) Math.max(0, dest.x - src.x);
        int rowLen = bpp * (int) Math.min(intersection.width, tileSizeX[regionIndex]);

        int outputRow = intersection.y - y;
        int outputCol = intersection.x - x;
        int outputOffset = outputRow * outputRowLen + outputCol * bpp;
        for (int c=0; c<getRGBChannelCount(); c++) {
          int srcChannelOffset = c * (tile.length / getRGBChannelCount());
          int destChannelOffset = c * w * h * bpp;
          for (int copyRow=0; copyRow<intersection.height; copyRow++) {
            int realRow = copyRow + intersection.y - src.y;
            int inputOffset = bpp * (realRow * tileSizeX[regionIndex] + intersectionX);
            System.arraycopy(tile, srcChannelOffset + inputOffset,
              buf, destChannelOffset + outputOffset + copyRow*outputRowLen, rowLen);
          }
        }
      }
    }
    catch (SQLException e) {
      LOGGER.warn("Failed to query tiles", e);
    }
    finally {
      try {
        conn.close();
      }
      catch (SQLException e) {
        LOGGER.warn("Failed to close connection", e);
      }
    }

    return buf;
  }

  // -- Internal FormatReader API methods --

  /* @see loci.formats.FormatReader#initFile(String) */
  @Override
  protected void initFile(String id) throws FormatException, IOException {
    super.initFile(id);

    findDBFiles();

    core.clear();

    JSONObject[] regionMetadata = new JSONObject[pixelsFiles.size()];
    fullResolutionCoreIndex = new int[pixelsFiles.size()];
    tileSizeX = new int[pixelsFiles.size()];
    tileSizeY = new int[pixelsFiles.size()];
    overlapX = new int[pixelsFiles.size()];
    overlapY = new int[pixelsFiles.size()];
    scaleFactor = new int[pixelsFiles.size()];
    tileRangeX = new int[2][pixelsFiles.size()];
    tileRangeY = new int[2][pixelsFiles.size()];
    zSteps = new Integer[pixelsFiles.size()][];

    Arrays.fill(tileRangeX[0], Integer.MAX_VALUE);
    Arrays.fill(tileRangeY[0], Integer.MAX_VALUE);

    for (int i=0; i<pixelsFiles.size(); i++) {
      String file = pixelsFiles.get(i);
      CoreMetadata m = new CoreMetadata();
      m.pixelType = FormatTools.UINT8;

      Connection conn = openConnection(file);
      try {
        PreparedStatement regionQuery = conn.prepareStatement(
          "SELECT id, data FROM region"
        );
        ResultSet regionData = regionQuery.executeQuery();
        // expect a single row returned
        if (regionData.next()) {
          String json = regionData.getString(2);
          LOGGER.trace("{}", json);

          try {
            // expect trailing whitespace and line breaks
            // in the value for "AcquisitionSettings",
            // which will prevent parsing
            json = json.trim();
            json = json.replaceAll("\r\n", "_");
            regionMetadata[i] = new JSONObject(json);
          }
          catch (JSONException je) {
            throw new FormatException("Could not read metadata for region in " + file, je);
          }

          // "Rows" and "Columns" in the JSON metadata here
          // reflect the canvas size, not the area (FOVs) actually acquired

          tileSizeX[i] = regionMetadata[i].getInt("ImageWidth");
          tileSizeY[i] = regionMetadata[i].getInt("ImageHeight");
          overlapX[i] = regionMetadata[i].getInt("OverlapWidth");
          overlapY[i] = regionMetadata[i].getInt("OverlapHeight");
          scaleFactor[i] = regionMetadata[i].getInt("CacheStep");
        }

        PreparedStatement fovQuery = conn.prepareStatement(
          "SELECT row, column FROM fovs"
        );
        ResultSet fovs = fovQuery.executeQuery();
        while (fovs.next()) {
          int row = fovs.getInt(1);
          int col = fovs.getInt(2);
          tileRangeY[0][i] = (int) Math.min(tileRangeY[0][i], row);
          tileRangeY[1][i] = (int) Math.max(tileRangeY[1][i], row);
          tileRangeX[0][i] = (int) Math.min(tileRangeX[0][i], col);
          tileRangeX[1][i] = (int) Math.max(tileRangeX[1][i], col);
        }
        // TODO: no overlap handling yet
        m.sizeX = (tileRangeX[1][i] - tileRangeX[0][i] + 1) * tileSizeX[i];
        m.sizeY = (tileRangeY[1][i] - tileRangeY[0][i] + 1) * tileSizeY[i];

        PreparedStatement maxLevelQuery = conn.prepareStatement(
          "SELECT level FROM images ORDER BY level DESC"
        );
        ResultSet maxLevel = maxLevelQuery.executeQuery();
        if (maxLevel.next()) {
          int resolutionCount = maxLevel.getInt(1);
          m.resolutionCount = resolutionCount + 1;
        }

        PreparedStatement channelQuery = conn.prepareStatement(
          "SELECT id, name, color, save_16bit, excitation_wavelength, emission_wavelength FROM channels ORDER BY id"
        );
        ResultSet channels = channelQuery.executeQuery();
        while (channels.next()) {
          m.sizeC++;

          // TODO: save channel name, color, wavelengths

          boolean save16 = channels.getBoolean(4);
          if (save16) {
            m.pixelType = FormatTools.UINT16;
          }
        }

        PreparedStatement zQuery = conn.prepareStatement(
          "SELECT DISTINCT is_zstack,z_position FROM images WHERE level=0 ORDER BY is_zstack,z_position"
        );
        ResultSet zs = zQuery.executeQuery();
        ArrayList<Integer> tmpZ = new ArrayList<Integer>();
        while (zs.next()) {
          boolean isZ = zs.getBoolean(1);
          int zPos = zs.getInt(2);

          tmpZ.add(zPos);
        }
        zSteps[i] = tmpZ.toArray(new Integer[tmpZ.size()]);
      }
      catch (SQLException e) {
        LOGGER.warn("Failed to initialize", e);
      }
      finally {
        try {
          conn.close();
        }
        catch (SQLException e) {
          LOGGER.warn("Failed to close connection", e);
        }
      }

      m.sizeZ = zSteps[i].length;
      m.sizeT = 1;
      m.imageCount = m.sizeZ * m.sizeC * m.sizeT;

      // TODO: bad assumption in general?
      if (m.sizeC == 1 && m.pixelType == FormatTools.UINT8) {
        m.sizeC = 3;
        m.rgb = true;
      }

      m.dimensionOrder = "XYCZT";

      core.add(m);
      fullResolutionCoreIndex[i] = core.size() - 1;
      for (int r=1; r<m.resolutionCount; r++) {
        CoreMetadata res = new CoreMetadata(m);
        int scale = (int) Math.pow(scaleFactor[i], r);
        res.sizeX /= scale;
        res.sizeY /= scale;
        res.resolutionCount = 1;
        core.add(res);
      }
    }

    // TODO: parse ROIs from .aqproj

    MetadataStore store = makeFilterMetadata();
    MetadataTools.populatePixels(store, this);

    String instrument = MetadataTools.createLSID("Instrument", 0);
    store.setInstrumentID(instrument, 0);

    for (int i=0; i<regionMetadata.length; i++) {
      int imageIndex = hasFlattenedResolutions() ? fullResolutionCoreIndex[i] : i;

      String objectiveID = MetadataTools.createLSID("Objective", 0, i);
      store.setObjectiveID(objectiveID, 0, i);

      Double lensNA = regionMetadata[i].getDouble("ObjectiveLensNA");
      store.setObjectiveLensNA(lensNA, 0, i);

      String immersion = regionMetadata[i].getString("ObjectiveImmersion");
      store.setObjectiveImmersion(getImmersion(immersion), 0, i);

      Double magnification = regionMetadata[i].getDouble("ObjectiveNominalMagnification");
      store.setObjectiveNominalMagnification(magnification, 0, i);

      String objectiveName = regionMetadata[i].getString("ObjectiveName");
      store.setObjectiveModel(objectiveName, 0, i);

      store.setImageName(regionMetadata[i].getString("Name"), imageIndex);
      store.setObjectiveSettingsID(objectiveID, imageIndex);

      Double physicalX = regionMetadata[i].getDouble("PhysicalSizeX");
      Double physicalY = regionMetadata[i].getDouble("PhysicalSizeY");

      store.setPixelsPhysicalSizeX(FormatTools.getPhysicalSizeX(physicalX), imageIndex);
      store.setPixelsPhysicalSizeY(FormatTools.getPhysicalSizeY(physicalY), imageIndex);
    }
  }

  // -- Helper methods --

  private void findDBFiles() {
    Location dir = new Location(getCurrentFile()).getAbsoluteFile().getParentFile();
    String[] list = dir.list(true);
    Arrays.sort(list);
    for (String f : list) {
      Location slide = new Location(dir, f);
      if (f.startsWith("Slide ") && slide.isDirectory()) {
        String[] dbFiles = slide.list(true);
        Arrays.sort(dbFiles);
        for (String db : dbFiles) {
          if (checkSuffix(db, "tfcyto")) {
            pixelsFiles.add(new Location(slide, db).getAbsolutePath());
          }
        }
      }
    }
  }

  private int getLevel() {
    if (hasFlattenedResolutions()) {
      int fullResIndex = fullResolutionCoreIndex[getRegionIndex()];
      return getCoreIndex() - fullResIndex;
    }
    return getResolution();
  }

  private int getRegionIndex() {
    if (hasFlattenedResolutions()) {
      int index = getCoreIndex();
      for (int i=fullResolutionCoreIndex.length-1; i>=0; i--) {
        if (fullResolutionCoreIndex[i] < index) {
          return i;
        }
      }
    }
    return getSeries();
  }

  private String getDBFile() {
    return pixelsFiles.get(getRegionIndex());
  }

  private Codec getCodec(int compression) throws UnsupportedCompressionException {
    switch (compression) {
      case 0:
      case 1:
        return new PassthroughCodec();
      case 6:
        return new JPEGXRCodec();
      case 7:
        return new JPEGCodec();
      default:
        throw new UnsupportedCompressionException("Unsupported compression: " + compression);
    }
  }

  private Connection openConnection(String file) throws IOException {
    Connection conn = null;
    try {
      // see https://github.com/xerial/sqlite-jdbc/issues/247
      SQLiteConfig config = new SQLiteConfig();
      config.setReadOnly(true);
      conn = config.createConnection("jdbc:sqlite:" +
        new Location(file).getAbsolutePath());
    }
    catch (SQLException e) {
      LOGGER.warn("Could not read from database");
      throw new IOException(e);
    }
    return conn;
  }

}
