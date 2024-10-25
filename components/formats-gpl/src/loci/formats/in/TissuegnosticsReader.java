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
import ome.xml.model.primitives.Color;
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

  // per specification, wavelengths outside this range should be ignored
  private static final int WAVE_MIN = 300;
  private static final int WAVE_MAX = 800;

  private List<String> pixelsFiles = new ArrayList<String>();
  private List<ScanRegion> regions = new ArrayList<ScanRegion>();

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
      regions.clear();
    }
  }

  /* @see loci.formats.IFormatReader#getOptimalTileWidth() */
  @Override
  public int getOptimalTileWidth() {
    FormatTools.assertId(currentId, true, 1);
    try {
      ScanRegion region = getRegion();
      return region.tileSizeX / (int) Math.pow(region.scaleFactor, getLevel());
    }
    catch (FormatException e) {
      LOGGER.warn("Could not get optimal tile width", e);
    }
    return super.getOptimalTileWidth();
  }

  /* @see loci.formats.IFormatReader#getOptimalTileHeight() */
  @Override
  public int getOptimalTileHeight() {
    FormatTools.assertId(currentId, true, 1);
    try {
      ScanRegion region = getRegion();
      return region.tileSizeY / (int) Math.pow(region.scaleFactor, getLevel());
    }
    catch (FormatException e) {
      LOGGER.warn("Could not get optimal tile height", e);
    }
    return super.getOptimalTileHeight();
  }

  /* @see loci.formats.IFormatReader#getSeriesUsedFiles(boolean) */
  @Override
  public String[] getSeriesUsedFiles(boolean noPixels) {
    FormatTools.assertId(currentId, true, 1);

    ArrayList<String> files = new ArrayList<String>();
    files.add(getCurrentFile());
    if (!noPixels) {
      try {
        files.add(getRegion().file);
      }
      catch (FormatException e) {
        LOGGER.warn("Could not get file", e);
      }
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

    int[] zct = getZCTCoords(no);
    ScanRegion region = getRegion(zct[2]);
    int bpp = FormatTools.getBytesPerPixel(getPixelType());
    int level = getLevel();
    int scale = (int) Math.pow(region.scaleFactor, level);

    int scaledOverlapX = region.overlapX / scale;
    int scaledOverlapY = region.overlapY / scale;

    Region dest = new Region(x, y, w, h);
    Connection conn = openConnection(region.file);
    try {
      PreparedStatement tiles = conn.prepareStatement(
        "SELECT data, compression, row, column FROM images WHERE region=? AND level=? AND " +
        "row>=? AND row<=? AND column>=? AND column<=? AND channel=? AND is_zstack=? AND z_position=? ORDER BY row,column"
      );

      int rowOffset = (int) Math.floor((double) region.tileRangeY[0] / scale);
      int colOffset = (int) Math.floor((double) region.tileRangeX[0] / scale);

      int startRow = (int) Math.floor((double) y / region.tileSizeY) + rowOffset;
      int startCol = (int) Math.floor((double) x / region.tileSizeX) + colOffset;
      int endRow = (int) Math.ceil((double) (y + h) / region.tileSizeY) + rowOffset;
      int endCol = (int) Math.ceil((double) (x + w) / region.tileSizeX) + colOffset;

      tiles.setInt(1, region.id);
      tiles.setInt(2, level);
      tiles.setInt(3, startRow);
      tiles.setInt(4, endRow);
      tiles.setInt(5, startCol);
      tiles.setInt(6, endCol);

      tiles.setInt(7, zct[1] + 1);
      tiles.setInt(8, zct[0] > 0 ? 1 : 0);
      tiles.setInt(9, region.zSteps[zct[0]]);

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
        options.width = region.tileSizeX;
        options.height = region.tileSizeY;

        Codec codec = getCodec(compression);
        byte[] tile = codec.decompress(data, options);

        int relativeColumn = column - (int) Math.floor((double) region.tileRangeX[0] / scale);
        int relativeRow = row - (int) Math.floor((double) region.tileRangeY[0] / scale);

        int pixelColumn = relativeColumn * options.width;
        int pixelRow = relativeRow * options.height;

        // overlap handling
        byte[][] fovs = null;
        Region[] fovPositions = new Region[scale * scale];
        if (level == 0) {
          Region fov = region.fovs.get(row + "-" + column);
          pixelRow -= fov.y;
          pixelColumn -= fov.x;

          pixelRow -= (relativeRow * region.overlapY);
          pixelColumn -= (relativeColumn * region.overlapX);

          fovPositions[0] = new Region(pixelColumn, pixelRow, options.width, options.height);
          fovs = new byte[1][];
          fovs[0] = tile;
        }
        else {
          fovs = splitFOVs(region, tile, scale);

          for (int f=0; f<fovs.length; f++) {
            int fovRow = row*scale + (f / scale);
            int fovColumn = column*scale + (f % scale);
            Region fov = region.fovs.get(fovRow + "-" + fovColumn);
            if (fov != null) {
              int xx = (fovColumn * options.width / scale) - (fov.x / scale) - (fovColumn * scaledOverlapX);
              int yy = (fovRow * options.height / scale) - (fov.y / scale) - (fovRow * scaledOverlapY);
              fovPositions[f] = new Region(xx, yy, options.width / scale, options.height / scale);
            }
          }
        }

        for (int f=0; f<fovs.length; f++) {
          if (fovPositions[f] == null) {
            continue;
          }
          Region intersection = fovPositions[f].intersection(dest);

          if (intersection.width > 0 && intersection.height > 0) {
            int outputRowLen = w * bpp;
            int intersectionX = (int) Math.max(0, dest.x - fovPositions[f].x);
            int rowLen = bpp * (int) Math.min(intersection.width, fovPositions[f].width);

            int outputRow = intersection.y - y;
            int outputCol = intersection.x - x;
            int outputOffset = outputRow * outputRowLen + outputCol * bpp;
            for (int c=0; c<getRGBChannelCount(); c++) {
              int srcChannelOffset = c * (fovs[f].length / getRGBChannelCount());
              int destChannelOffset = c * w * h * bpp;
              for (int copyRow=0; copyRow<intersection.height; copyRow++) {
                int realRow = copyRow + intersection.y - fovPositions[f].y;
                int inputOffset = bpp * (realRow * (region.tileSizeX / scale) + intersectionX);
                System.arraycopy(fovs[f], srcChannelOffset + inputOffset,
                  buf, destChannelOffset + outputOffset + copyRow*outputRowLen, rowLen);
              }
            }
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

    for (int i=0; i<pixelsFiles.size(); i++) {
      String file = pixelsFiles.get(i);
      CoreMetadata m = new CoreMetadata();
      m.pixelType = FormatTools.UINT8;

      int startRegionIndex = regions.size();

      Connection conn = openConnection(file);
      try {
        PreparedStatement regionQuery = conn.prepareStatement(
          "SELECT id, data, is_timelapse FROM region ORDER BY id"
        );
        ResultSet regionData = regionQuery.executeQuery();

        // expect one row per timepoint
        while (regionData.next()) {
          int regionID = regionData.getInt(1);
          String json = regionData.getString(2);
          boolean isTimelapse = regionData.getBoolean(3);

          LOGGER.trace("{}", json);

          ScanRegion region = new ScanRegion();
          region.id = regionID;
          region.file = file;

          try {
            // expect trailing whitespace and line breaks
            // in the value for "AcquisitionSettings",
            // which will prevent parsing
            json = json.trim();
            json = json.replaceAll("\r\n", "_");
            region.regionMetadata = new JSONObject(json);
          }
          catch (JSONException je) {
            throw new FormatException("Could not read metadata for region in " + file, je);
          }

          region.parseJSON();
          region.timepoint = regions.size() - startRegionIndex;

          // TODO: this means TMA regions are parsed, but not recorded separately
          // that might be OK, but needs to be double-checked
          if (isTimelapse || regions.size() == startRegionIndex) {
            regions.add(region);
          }
        }
        int timepoints = regions.size() - startRegionIndex;

        for (int regionIndex=startRegionIndex; regionIndex<regions.size(); regionIndex++) {
          ScanRegion currentRegion = regions.get(regionIndex);

          PreparedStatement fovQuery = conn.prepareStatement(
            "SELECT row, column, stitch_rectangle_x, stitch_rectangle_y, stitch_rectangle_w, stitch_rectangle_h FROM fovs WHERE region_id=?"
          );
          fovQuery.setInt(1, regions.get(regionIndex).id);

          ResultSet fovs = fovQuery.executeQuery();
          while (fovs.next()) {
            int row = fovs.getInt(1);
            int col = fovs.getInt(2);
            double x = fovs.getDouble(3);
            double y = fovs.getDouble(4);
            double w = fovs.getDouble(5);
            double h = fovs.getDouble(6);
            currentRegion.tileRangeY[0] = (int) Math.min(currentRegion.tileRangeY[0], row);
            currentRegion.tileRangeY[1] = (int) Math.max(currentRegion.tileRangeY[1], row);
            currentRegion.tileRangeX[0] = (int) Math.min(currentRegion.tileRangeX[0], col);
            currentRegion.tileRangeX[1] = (int) Math.max(currentRegion.tileRangeX[1], col);

            Region fov = new Region((int) x, (int) y, (int) w, (int) h);
            currentRegion.fovs.put(row + "-" + col, fov);
          }
          int xTiles = currentRegion.tileRangeX[1] - currentRegion.tileRangeX[0] + 1;
          int yTiles = currentRegion.tileRangeY[1] - currentRegion.tileRangeY[0] + 1;
          m.sizeX = xTiles * (currentRegion.tileSizeX - currentRegion.overlapX);
          m.sizeY = yTiles * (currentRegion.tileSizeY - currentRegion.overlapY);

          PreparedStatement maxLevelQuery = conn.prepareStatement(
            "SELECT level FROM images WHERE region=? ORDER BY level DESC"
          );
          maxLevelQuery.setInt(1, currentRegion.id);
          ResultSet maxLevel = maxLevelQuery.executeQuery();
          if (maxLevel.next()) {
            int resolutionCount = maxLevel.getInt(1);
            m.resolutionCount = resolutionCount + 1;
          }

          PreparedStatement zQuery = conn.prepareStatement(
            "SELECT DISTINCT is_zstack,z_position FROM images WHERE region=? ORDER BY is_zstack,z_position"
          );
          zQuery.setInt(1, currentRegion.id);

          ResultSet zs = zQuery.executeQuery();
          ArrayList<Integer> tmpZ = new ArrayList<Integer>();
          while (zs.next()) {
            boolean isZ = zs.getBoolean(1);
            int zPos = zs.getInt(2);

            tmpZ.add(zPos);
          }
          currentRegion.zSteps = tmpZ.toArray(new Integer[tmpZ.size()]);
          currentRegion.fullResolutionCoreIndex = core.size();
        }

        PreparedStatement channelQuery = conn.prepareStatement(
          "SELECT DISTINCT id, name, color, save_16bit, excitation_wavelength, emission_wavelength FROM channels ORDER BY id"
        );
        ResultSet channels = channelQuery.executeQuery();
        while (channels.next()) {
          m.sizeC++;

          Channel ch = new Channel();
          ch.id = channels.getInt(1);
          ch.name = channels.getString(2);
          ch.color = channels.getInt(3);

          boolean save16 = channels.getBoolean(4);
          if (save16) {
            m.pixelType = FormatTools.UINT16;
          }

          ch.exWave = channels.getInt(5);
          ch.emWave = channels.getInt(6);

          regions.get(startRegionIndex).channels.add(ch);
        }
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

      m.sizeZ = regions.get(startRegionIndex).zSteps.length;
      m.sizeT = regions.size() - startRegionIndex;
      m.imageCount = m.sizeZ * m.sizeC * m.sizeT;

      // TODO: bad assumption in general?
      if (m.sizeC == 1 && m.pixelType == FormatTools.UINT8) {
        m.sizeC = 3;
        m.rgb = true;
      }

      m.dimensionOrder = "XYCZT";
      m.littleEndian = true;

      core.add(m);
      for (int r=1; r<m.resolutionCount; r++) {
        CoreMetadata res = new CoreMetadata(m);
        int scale = (int) Math.pow(regions.get(startRegionIndex).scaleFactor, r);
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

    for (int i=0, index=0; i<regions.size(); index++) {
      ScanRegion region = regions.get(i);
      int imageIndex = hasFlattenedResolutions() ? region.fullResolutionCoreIndex : index;

      String objectiveID = MetadataTools.createLSID("Objective", 0, index);
      store.setObjectiveID(objectiveID, 0, index);

      Double lensNA = region.regionMetadata.getDouble("ObjectiveLensNA");
      store.setObjectiveLensNA(lensNA, 0, index);

      String immersion = region.regionMetadata.getString("ObjectiveImmersion");
      store.setObjectiveImmersion(getImmersion(immersion), 0, index);

      Double magnification = region.regionMetadata.getDouble("ObjectiveNominalMagnification");
      store.setObjectiveNominalMagnification(magnification, 0, index);

      String objectiveName = region.regionMetadata.getString("ObjectiveName");
      store.setObjectiveModel(objectiveName, 0, index);

      store.setImageName(region.regionMetadata.getString("Name"), imageIndex);
      store.setObjectiveSettingsID(objectiveID, imageIndex);

      Double physicalX = region.regionMetadata.getDouble("PhysicalSizeX");
      Double physicalY = region.regionMetadata.getDouble("PhysicalSizeY");

      store.setPixelsPhysicalSizeX(FormatTools.getPhysicalSizeX(physicalX), imageIndex);
      store.setPixelsPhysicalSizeY(FormatTools.getPhysicalSizeY(physicalY), imageIndex);

      for (int c=0; c<region.channels.size(); c++) {
        Channel ch = region.channels.get(c);
        store.setChannelName(ch.name, imageIndex, c);

        if (ch.emWave >= WAVE_MIN && ch.emWave <= WAVE_MAX) {
          store.setChannelEmissionWavelength(FormatTools.getWavelength((double) ch.emWave), imageIndex, c);
        }
        if (ch.exWave >= WAVE_MIN && ch.exWave <= WAVE_MAX) {
          store.setChannelExcitationWavelength(FormatTools.getWavelength((double) ch.exWave), imageIndex, c);
        }
        store.setChannelColor(ch.getColor(), imageIndex, c);
      }

      i += core.get(imageIndex).sizeT;
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

  private ScanRegion getRegion() throws FormatException {
    return getRegion(0);
  }

  private ScanRegion getRegion(int t) throws FormatException {
    int index = getCoreIndex();
    for (int i=regions.size()-1; i>=0; i--) {
      ScanRegion r = regions.get(i);
      if (r.timepoint == t && r.fullResolutionCoreIndex <= index) {
        return r;
      }
    }
    throw new FormatException("Could not find ScanRegion (core index " + index + ", t=" + t + ")");
  }

  private int getLevel() throws FormatException {
    if (hasFlattenedResolutions()) {
      ScanRegion region = getRegion();
      return getCoreIndex() - region.fullResolutionCoreIndex;
    }
    return getResolution();
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

  /**
   * The largest resolution stores one field of view (FOV) per tile.
   * Each sub-resolution has the same tile size as the largest resolution,
   * so the number of FOVs stored in one tile increases as the resolutions
   * get smaller.
   *
   * For example, with a tile size of 2048x2048 and a downsample factor of 4,
   * resolution 1 will have 16 (4x4) FOVs per tile, with each FOV being 512x512.
   *
   * However, none of the FOV positions or overlaps are taken into account.
   * So the first step in assembling a sub-resolution is to split each tile
   * into its component FOVs. Then each FOV can be repositioned according to
   * its associated stitching rectangle and overlap.
   *
   * This method only splits the given tile into component FOVs.
   * This is not necessary for the largest resolution, only sub-resolutions.
   */
  private byte[][] splitFOVs(ScanRegion region, byte[] tile, int scale) {
    byte[][] fovs = new byte[scale * scale][];
    int bpp = FormatTools.getBytesPerPixel(getPixelType());
    int channels = getRGBChannelCount();

    int srcWidth = region.tileSizeX * bpp;
    int destWidth = (region.tileSizeX / scale) * bpp;
    int srcHeight = region.tileSizeY;
    int destHeight = region.tileSizeY / scale;

    for (int fov=0; fov<fovs.length; fov++) {
      int fovRow = fov / scale;
      int fovCol = fov % scale;

      fovs[fov] = new byte[destWidth * destHeight * bpp * channels];

      for (int c=0; c<channels; c++) {
        int srcChannelOffset = c * srcWidth * srcHeight;
        int destChannelOffset = c * destWidth * destHeight;

        for (int row=0; row<destHeight; row++) {
          int srcOffset = srcChannelOffset + (((fovRow * destHeight) + row) * srcWidth) + (fovCol * destWidth);
          int destOffset = destChannelOffset + (row * destWidth);

          System.arraycopy(tile, srcOffset, fovs[fov], destOffset, destWidth);
        }
      }
    }

    return fovs;
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

  class ScanRegion {
    public String file;
    public JSONObject regionMetadata;
    public int id;
    public int fullResolutionCoreIndex;
    public int tileSizeX;
    public int tileSizeY;
    public int overlapX;
    public int overlapY;
    public int scaleFactor;
    public int[] tileRangeX = new int[] {Integer.MAX_VALUE, 0};
    public int[] tileRangeY = new int[] {Integer.MAX_VALUE, 0};
    public Integer[] zSteps;
    public int timepoint;
    public List<Channel> channels = new ArrayList<Channel>();
    public HashMap<String, Region> fovs = new HashMap<String, Region>();

    public void parseJSON() {
      // "Rows" and "Columns" in the JSON metadata here
      // reflect the canvas size, not the area (FOVs) actually acquired
      tileSizeX = regionMetadata.getInt("ImageWidth");
      tileSizeY = regionMetadata.getInt("ImageHeight");
      overlapX = regionMetadata.getInt("OverlapWidth");
      overlapY = regionMetadata.getInt("OverlapHeight");
      scaleFactor = regionMetadata.getInt("CacheStep");
    }

  }

  class Channel {
    public int id;
    public String name;
    public int color;
    public int exWave;
    public int emWave;

    /** Color returned by DB is ARGB, OME model is RGBA. */
    public Color getColor() {
      int alpha = (color >> 24) & 0xff;
      int red = (color >> 16) & 0xff;
      int green = (color >> 8) & 0xff;
      int blue = color & 0xff;
      return new Color(red, green, blue, alpha);
    }
  }

}
