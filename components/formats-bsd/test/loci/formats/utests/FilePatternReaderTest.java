/*
 * #%L
 * BSD implementations of Bio-Formats readers and writers
 * %%
 * Copyright (C) 2005 - 2017 Open Microscopy Environment:
 *   - Board of Regents of the University of Wisconsin-Madison
 *   - Glencoe Software, Inc.
 *   - University of Dundee
 * %%
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * #L%
 */

package loci.formats.utests;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertFalse;
import static org.testng.Assert.assertTrue;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;


import loci.formats.IFormatReader;
import loci.formats.in.FilePatternReader;
import loci.formats.Memoizer;
import loci.formats.ome.OMEXMLMetadata;
import loci.formats.meta.MetadataRetrieve;
import loci.formats.services.OMEXMLService;
import loci.common.services.ServiceFactory;

import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;


public class FilePatternReaderTest {

  private Path path;
  private String id;
  private FilePatternReader reader = new FilePatternReader();
  private Memoizer memoizer;
  private OMEXMLService service;
  private MetadataRetrieve m;


  @BeforeMethod
  public void setUp() throws Exception {
    path = Files.createTempFile("test", ".pattern");
    id = path.toFile().getAbsolutePath();
    System.out.println(id);
    path.toFile().deleteOnExit();
    
    byte[] buf = "test_t<1-10>_c<0-3>_z<1-5>.fake".getBytes();
    Files.write(path, buf);
  }

  @AfterMethod
  public void tearDown() throws Exception {
    reader.close();
    memoizer.close();
  }

  public void assertDimensions() throws Exception {
    assertEquals(reader.getSizeX(), 512);
    assertEquals(reader.getSizeX(), 512);
    assertEquals(reader.getSizeT(), 10);
    assertEquals(reader.getSizeC(), 4);
    assertEquals(reader.getSizeZ(), 5);
  }

  @Test
  public void testSimple() throws Exception {
    reader.setId(id);
    assertDimensions();
    reader.close();
    reader.setId(id);
    assertDimensions();
  }

  @Test
  public void testMemoizer() throws Exception {
    memoizer = new Memoizer(reader);
    memoizer.setId(id);
    assertDimensions();
    memoizer.close();
    memoizer.setId(id);
    assertDimensions();
  }
}
