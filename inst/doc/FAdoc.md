---
author:
- |
    Alex Deckmyn\
    Royal Meteorological Institute of Belgium,\
    Ringlaan 3, B-1180 Brussels, Belgium,\
    e-mail: alex.deckmyn@meteo.be
title: |
    File formats in the ALADIN universe\
    The internal structure of LFI, FA, Surfex and MF files 
---

**THIS IS WORK IN PROGRESS!!!**

# Introduction

This document grew from re-coding FA file access in R. When working on
this code, one problem was the limited documentation on the actual format
of the files.

In parallel to this document, an implementation in R and C was written
for the various file formats. Therefore, this code may serve as a
further documentation of the internal structure of the files. However,
it does not implement all possible flavours and exceptions.\

**AKNOWLEDGEMENTS:**\

In preparing this document and the code implementing the various
formats, I have found inspiration in the following documents & code:

-   LFI and FA schematic descriptions by Jean-Daniël Gril
    (Météo France).

-   FA library documentation (cy24) by J. Clochard, R. El Katib & D.
    Paradis (Météo France).

-   Basic LFI implementation in the GL package by Ulf Andrae (Sweden)

-   The standard FA implementation in the ALADIN export version.

-   Some documentation from Meso-NH for the SURFEX style LFI files.

There may be more documentation somewhere, but it does not appear to be
readily available to the ALADIN community.\

**WARNING & INVITATION:**\

-   This document focusses on formats as they are currently in use.
    Therefore, some early (pre-1990) cases may not be included.

-   Not all aspects of the LFI information headers are described. Most
    of them are for “bookkeeping” and are not vital to reading the
    data files. My basic code does not use them.

-   Therefore, if you use different tools than the standard library, you
    may create files that violate some specification that I don’t
    know about. If you’re only reading data, that shouldn’t worry.

-   This document is still evolving. Some special cases may not be
    described yet. Please feel free to contribute!

#LFI: the basic structure

All ALADIN, Arpège and SURFEX files are based on this format. It is very
generic and may in principle contain any kind of data whatsoever. The
graphical representation by J-D Gril gives a very clear and compact
overview of the full structure. This section merely adds some details &
personal considerations.

##Conventions

-   All numbers (real and integer) are supposed to be BIG-ENDIAN.

-   LFI is based on a word length of 8 octets. Specifically, this means
    that in the 3 “header” sections all integers are 64-bit (signed? I
    guess so. Certainly in the data sections.). The actual encoded data
    may have any format.

-   Addresses and sector lengths are also given in words (64 bits),
    not octets.

-   *In `lficom0.h` the word length is fixed to 8 except on
    HPPA platform. I don’t know what this would have as an effect on
    the files.*

##The 4 sections structure

An LFI file consists of 4 sections. The first three are “headers” while
the fourth contains the actual binary data. The three header sections
have exactly the same length. In the original Fortran code, the file is
read in “DIRECT ACCESS” mode, requiring a fixed length for all records.
The first 3 sections are exactly 1 record, while the data section can be
much longer (but always a multiple of the record size).

1.  **Information section.** Usually this section contains 22 words
    of meta-data. The first word contains the *length of a physical
    records* (sector length SL), expressed in words of 8 octets. The
    length of sections 1 to 3 in octets is exactly 8\*SL. A typical
    value is `C00h`=3072. Yes, this means that section 1 has a lot of
    open space! It should be filled with 0, but you could write a poem
    or even a short novella (24000 ASCII characters) there and no-one
    would ever know. Maybe some people do this already?

2.  **Article name section.** A list of 16-character (=2 words) names
    for the data articles. The length (16) is given in the second word
    of section 1, but if it isn’t 16 we’re in trouble. Shorter names
    should be padded with ASCII space, not 0 bits! The names are not \\0
    terminated like in C.

3.  **Article length and location section.** Contains a list of length
    and position in the file (in words) for every article.

4.  **Data section.** Variable length. Every article begins at a new
    word, so there may be some trailing bits at the end of an article.
    The data section may consists of many (fixed-size) records and
    articles are not restricted in any way by the boundaries of
    the records.

Note that in the sections 2 and 3, you need 16 octets per data field.
This also means that the record length limits the number of fields that
can be encoded in the name sector. To get around that restriction an
expansion mechanism was introduced (see further). This is very different
from e.g. GRIB, where new data can just be attached to the file at any
time.

##Section 1: information

This section only uses 22 words out of (typically) almost 20KB. And even
then, only one or two are really important for decoding. The meaning of
the entries is given `lficom0.h`. I’ll translate “Article physique” as a
*data sector* (in fact the fortran name would be a record!) and “Article
Logique” as a *data record* (either a data field or part of the frame).

-   S1[1] This is the “length of a sector” (sector length SL),
    expressed in words of 8 octets (64 bits). In fortran, this would be
    called the record length (except that it may not be expressed in 64
    bit words). The first 3 sections have exactly this length. The data
    section can be much longer, but always a multiple of this. Data
    fields don’t have to be in separate physical sectors. SL must be
    even and the number of data fields is restricted to SL/2, as
    explained above.

-   S1[2] Maximum length of article names (should always be 16).

-   S1[3] Flag for correct closure of file (should always be 0.

-   S1[4] Length of sector 1 (should always be 22).

-   S1[5] total number of data sectors in the file. The header
    sections account for 3.

-   S1[6] Number of data records in the file. Note that in an
    ALADIN/Arpege file, the first 7 records form the *frame* and
    concentrate all the meta data. In other files, e.g. Surfex, the
    number of “meta-data” fields is not fixed, and they can occur
    anywhere in the list. Also, some of the records counted here may in
    fact be “holes” (explained below).

-   S1[7] Length of shortest data record (in FA file: should be 1,
    because the 7th (frame) article has length 1.

-   S1[8] Length of longest data record.

-   S1[9] Total length of data.

-   S1[10] Number of in-place re-writes (=overwrite data with
    same length)

-   S1[11] Number of re-writes with shorter data.

-   S1[12] Number of re-writes with longer data (these will leave the
    original location as a hole.)

-   S1[13] Maximum number of data records per “sequence”
    (=sector length/2)

-   S1[14] Creation date (integer yyyymmdd).

-   S1[15] Creation hour (integer hhmmss).

-   S1[16] Date of last modification (integer yyyymmdd).

-   S1[17] Hour of last modification (integer hhmmss).

-   S1[18] Date of first modification.

-   S1[19] Hour of first modification.

-   S1[20] Number of pairs of index sectors. Usually 1. I don’t really
    know what this means. *It is 1 even if there is a second sector for
    file names*.

-   S1[21] Number of holes due to rewrites.

-   S1[22] Maximum number of data sectors used.

Many of these entries are more for bookkeeping (date of creation, date
of last modification etc.). You can read all information in the file
without ever referencing them. If you change something in the file,
though, you may want to do the bookkeeping.

The SL is chosen when creating the file. A typical value is 3072, but
this is not a rule. As far as I understand, this length should be long
enough because SL/2 is the maximal number of logical articles (fields)
that the file can describe (2 words per article name) without expanding
it. If you want to add fields to an existing file, it is easy to do this
until you fill the header sections completely. Then you must create new
name and address sectors (next paragraph).

If S[6]$>$S[13], there is a **second set of name & address
sectors**. To find this second sector with continuation of the field
names, the address is encoded **at the end of the header section**. In
reverse order (I think), the sector numbers of the next name sectors are
given as 8 octet integers. Note that this is not an offset, but an
address in units of S1[1]. So in C code, you would do something like

    fseek(fafile, sector_size * 8 - (next_sector_number - 1) * 8, SEEK_SET);
    // read 8 byte integer next_sector
    fseek(fafile, sector_size * 8 * (next_sector - 1),SEEK_SET);

##Section 2: article names

The names of the logical articles are written as ASCII characters. The
maximal length of a name (in bytes, not words this time!) is given by
S1[2]. Usually it should be 16. If it is different, I have no idea
what would happen. Luckily, I have never seen this.

If a name is shorter than 16 characters, it should be padded with space
characters (not “0” but the ASCII value for blank space, i.e. 32).
Names are not finalised with \\0 like a string in C. So you probably
shouldn’t attempt to read this section with something like
`readf(file,**string);`

The locations that have not yet been used, are not filled with 0 or
blanks, but with the 16 character word `**FIN D’INDEX** `. But as with
the first section, you can probably have some fun in this unused part of
the file.

##Section 3: length and offset of data records

For every data record you get 2 words (64 bit integers). The first gives
the length (in words) of the data, the second gives the location in the
file. Note that this second value, the offset, starts from 0 and is also
given in words. So offset=1000 means that the data starts from the
1000st word in the file. Depending on your code (C, fortran, R, python…)
this may either make you very happy or be the cause of some serious
headaches. Notice that the above means that the the length of a data
record is always given in words, so the number of bytes is a multiple of
8.

Anyway, in C something like `fseek(file,8*(offset-1))` should leave you
at the right spot to start reading the data. It’s much harder in
fortran. In fact, in my opinion fortran is not at all a suitable
language for accessing binary files.

##Section 4: data

There are no rules here! A data record may contain anything. ASCII text,
a GRIB file, a binary dump… It is completely up to the user (well OK, the
developer) to know how to read the data.

Data articles may extend over multiple physical sector (SL). And
conversely, there may be several shorter articles within 1 physical
sector. So I don’t know how far these SL-sized building blocks have any
significance in the data section. They probably are just a signature of
the original fortran code. For fortran routines, it is vital that every
LFI file has a total length that is exactly a multiple of the sector
length!

##Some considerations

-   Adding a new article to an existing LFI file is simple. You add the
    name, length and offset, you increase the article counter and add
    the data at the end. Et voilà! If you have to add a data sector at
    the end, you should also change the total file length and increase
    the sector counter.

-   Replacing an existing article is also simple if the new version has
    the same dimension. You just overwrite (OK, and then you update all
    the counters and dates in the header if you’re nice).

-   If the replacement article is shorter than the original, than you
    can still overwrite it and you just have to modify the article
    length entry in section 3. Yes, this leaves a hole in the file. But
    no worries, the code doesn’t have to know. And if you really like
    bit-counting, there even is a bookkeeping counter for this in
    section 1. You could re-align all the data, but that is hard work.

-   If the replacement is longer than the original (e.g. by changing the
    GRIB compression ratio), you have a bit more work. Either you
    re-align all the next articles (by which I mean: change the offset
    address to fit the new data). Or you add the new version at the end
    of the list and erase the first one in sections 2 and 3. That leaves
    a big hole, of course. These are counted in S1(21). And S1(5), the
    field counter, actually includes these holes. These holes are
    indicated simply by setting the article name to ’ ’ and may be
    re-used for a new field.

-   This means that you have to know exactly how long your encoded data
    is, before you can decide where it will be placed in the
    data section.

**NOTE:** You can in fact add more fields to a file than allowed by the
length of the name sector. It appears that expanding an LFI/FA file with
a second (third …) name sector is done in the following order:

1.  First write the data as usual, just after the previous data set.

2.  Notice that the name sector is full, so you need a new one.

3.  Fill the current data sector with zeroes and then start new name &
    address sectors, which have one entry.

4.  Write the sector number of the new name sector to the **end** of the
    header sector.

As a result, the first field of the next “extension” usually comes
*before the actual name sector*.

#LFI SURFEX files

In this section on MesoNH files, I will only discuss the files from
SURFEX that are relevant to ALADIN. Because of that, I will not discuss
compression (for which I have no documenation anyway) as it does not
appear to be relevant to SURFEX output.

Typically, SURFEX output files are called `AROMOUT_.00hh.lfi`, where
`hh` denotes the forecast leadtime. All 2D fields are in gridpoint
space. The basic decoding routines from GL were very helpful in
understanding the data structure.

##Meta-data

All information about grid dimension, projection etc. is coded into
different logical articles. Some of these contain only 1 single number
(real or integer). It is important to know whether they are long
integers of doubles. This is not indicated explicitely in the file
(unless maybe in the comment section described below). The structure of
the data articles is discussed in the next subsection.

The most important ones are:

-   <span>**IMAX,JMAX**</span>: These are the grid dimensions.

-   <span>**CURR%TDATA**,CURR%TIME</span>: The current date & time (i.e.
    forecast time + lead time)

##Data organisation

All LFI articles in an FM file are organised as follows:

1.  Word 1: grid type (0 for a scalar, 2 for a 2D grid point field,
    beyond that I don’t know)

2.  Word 2: length of the comment section (in words, i.e. 8 octets). By
    definition, every article should start with a 100 octet comment
    section, so this should always be 100.

3.  Word 3-102: comment section. Strangely, this comment is stored as
    `long integers` (64 bit) but every integer represents one
    ASCII character. So you should read them as long integers and then
    re-cast them as character data. Or maybe this is meant to allow
    non-ASCII symbols, too? I guess you could fit Kanji if you wanted.
    Anyway, be careful.

4.  Remainder: data. The length of this data is to be calculated as the
    length read from section 3 (see previous chapter) minus (2 +
    comment length). For a 2D field, this data is simply a list of
    8-octet doubles. Some fields are coded as integer (no problem there)
    or even as logical. In the latter case, you need only 1 bit per
    grid point. I haven’t figured out yet exactly how the data is
    organised (e.g. is the data padded with zeroes to fill a number
    of words?). But this data is not so important to us. **Note that
    data for an N x M grid is stored as a (N+2)x(M+2) matrix. So after
    getting the data in the right matrix shape, you should still peal
    off one layer of grid points. To add to the confusion, the SW point
    given in the metadata is in fact in this unphysical part, while the
    grid dimension is given for the physical grid. Not
    very consistent...**

##Example code in R

For example, this set of commands reads the typical fields necessary for
defining the projection:

      LFIproj <- trim(int2ascii(LFIread(lfi,"GRID_TYPE",type="integer")))

      if(LFIproj=="CONF PROJ"){
        nx <- LFIread(lfi,"IMAX","integer")
        ny <- LFIread(lfi,"IMAX","integer")
        SW <- c(LFIread(lfi,"LONOR","numeric"), LFIread(lfi,"LATOR","numeric"))
        Lat0 <- LFIread(lfi,"LAT0","numeric")
        Lon0 <- LFIread(lfi,"LON0","numeric")
    #  rpk <-  LFIread(lfi,"RPK","numeric")  ### not necessary
    #  beta <- LFIread(lfi,"BETA","numeric") ### assumed 0
        xhat <- LFIread(lfi,"XHAT","numeric") 
        yhat <- LFIread(lfi,"YHAT","numeric")
        dx <- xhat[2] - xhat[1]
        dy <- yhat[2] - yhat[1]
        SW <- SW + c(dx,dy)
        projection=list(proj="lcc",lon_0=Lon0,lat_1=Lat0,lat_2=Lat0,a=6371229.0,es=0.0)

Notice how some values are real and others are integer or even
characters. The field `GRID_TYPE` will usually have the value
`CONF PROJ` for Lambert conformal projection.

Also notice how we add (dx,dy) to SW. The reason is that, as far as I
can see, the SW co-ordinates given in the file are for the “unphysical”
extension of the grid.

#FA files

Alaro/Arome files are organised differently. They still follow the LFI
standard, though. So the only difference is in the logical articles (the
data).

##The frame

First of all, the geometry, date, encoding… of the data are gathered in 7
logical articles (whereas in FM almost every parameter was separate).
This is why you always have 7 **(since about cy40: 8!)** logical
articles more than the actual data fields. They should be the first
articles in the list, and they are:

1.  `CADRE_DIMENSIONS` 5 integers

2.  `CADRE_FRANKSCHMI` 4 reals

3.  `CADRE_REDPOINPOL` $8+2*(\mbox{NSMAX}+2)$ integers

4.  `CADRE_SINLATITUD` X reals

5.  `CADRE_FOCOHYBRID` $1+2*(1+NLEV)$ reals

6.  `Frame name`: the name is variable and indicates the name of
    this frame. Not very important in analysis. The data sector contains
    1 integer indicating the FA version, which should always be 1.

7.  `DATE-DES-DONNEES` 11 integers.

8.  `DATX-DES-DONNEES` 11 integers.

The first 5 articles define the **frame** (french: *cadre*), the
geometry (projection, hybrid levels, dimensions etc).

Looking at the routine `FAITOU`, it appears that the “frame name” field
**must** be the field following `CADRE_FOCOHYBRID`. That is the only way
to identify it, as the name is not fixed. Purely code-wise, for the rest
the frame fields could be anywhere and `FAITOU` would still work. But
they are *supposed* to be the first 5 fields, and the date is field 7.

The meaning of all the parameters is taken from the FA library
documentation. Note that several of the parameter names are derived from
their meaning in Arpège. For ALADIN files, the 7 (8) frame articles have
the same dimensions, but might sometimes mean something completely
different!

### CADRE\_DIMENSIONS

                              ALADIN                                           ARPEGE
  ------- ----------------------------------------------- ------------------------------------------------
   int 1      NSMAX=spectral truncation along Y axis                      NSMAX=truncation
   int 2        NDGL=number of points along Y axis         NDGL=number of latitude rows from pole to pole
   int 3        NDLON=number of points along X axis                            NDLON
   int 4         NFLEVG=number of vertical levels                              NFLEVG
   int 5   (-1)\*NMSMAX=spectral truncation along X axis      NSTTYP=kind of horizontal transformation

Notice that for ALADIN files, the last entry is negative. This serves as
a flag to differentiate ALADIN and ARPEGE files!

### CADRE\_FRANKSCHMI

In Arpege, this set of 4 reals describes the *stretching* of the
spectral data. In Aladin, only the first number is used to signify that
the domain is rotated. The fourth may be used to signal an academic case
(a torus with no extension zone), but I’ve never seen that either.

In modern ALADIN, the rotation angle BETA is only used for the mercator
projection (but I don’t know where...). In Lambert projection, the
rotation is given by LON0 (pole rotation is obsolete).

           ALADIN                                                 ARPEGE
  -------- ------------------------------------------------------ ----------------------------------------------
   real 1  BETA=angle of rotation                                 sine of the latitude of the pole of interest
   real 2  not used                                               cosine of longitude
   real 3  not used                                               sine of longitude
   real 4  PCODIL: 0 for real domain, 1 for virtual toric (???)   stretching factor

In most ALADIN cases, you may expect this section to contain only 0’s.

### CADRE\_REDPOINPOL

A list of 8+2\*(NSMAX+2) integers. The first 8 are important, because
they define the interpolation zone of the domain and some other
characteristics:

          ALADIN                                                  ARPEGE
  ------- ------------------------------------------------------- --------
   int 1  subtruncation                                           
   int 2  0:C+I gp,
          -1:C+I+E gp,
          1:C+I+E spectral for dynamical fields and gp for rest  
   int 3  NDLUX (usually 1)                                       
   int 4  NDLUN                                                   
   int 5  NDGUX (1)                                               
   int 6  NDGUN                                                   
   int 7  I zone in X axis (8)                                    
   int 8  I zone in Y axis (8)                                    

1.  Packed subtruncation

2.  domain representation: 0 (gridpoint C+I), -1 (gridpoint C+I+E), 1
    (spectral for dynamical fields, rest is gridopint on C+I+E)

3.  NDLUX (usually 1)

4.  NDLUN

5.  NDGUX

6.  NDGUN

7.  dimension of relaxation zone along X (typically 8)

8.  dimension of relaxation zone along Y (typically 8)

The rest of this vector is not used *according to the FA documentation*,
but in fact it holds values representing the number of spectral
components per line (according to the elliptic truncation):
$$cumsum(floor(nsmax*sqrt(abs(1-M^2/nmsmax^2))+1e-10+1)*4)$$ These **are
used in the XRD code** (contrary to the documentation!) when mixing the
compacted and non-compacted values, as explained later. Personally, I
found it easier to calculate these numbers whenever needed, rather than
having to get them from the frame. But Arpège may be different.

### CADRE\_SINLATITUD

In ARPEGE, this article gives a list of the sine of the latitude rows
for the northern hemisphere. It is a real vector with length at least
half the number of latitudes, and the length must be even.

In ALADIN files, however, it is a vector of (exactly) 18 reals
describing the projection and geometry of the domain. To be exact, since
about 2005 (“new eggx”) only 16 reals are used (and in a different
order!). So we must aditionally find out what version of ALADIN
formatting is used (you would think they could have changed the FA
version in the `Frame name` article, but that would of course be much
too simple). In eggx, this is the vector ZSINLA.

The distinction between “old style” and “new style” ALADIN files is made
by the first element. If it is zero or positive, we have an “old style”
file. For “new style” geometry, this element is -1 or -2. **This
explanation is taken from J-D Gril’s documentation. In the older FA code
documentation, this parameter is either 0 (no rotation) or 1 (the point
PLONR,PLATR is rotated onto the equator). I DO NOT KNOW WHETHER THIS IS
STILL VALID!!!** In any case, PLONR and PLATR are no longer available in
the “new style”. So it seems this +1 value is indeed obsolete.

Below I give a quick overview. More details can be found in J-D Gril’s
document, from which I got my information.

             **old**       **new**
  --------- ---------- ---------------
   real 1     NROTEQ       NROTEQ
   real 2     PLONR         ERPK
   real 3     PLATR         ELON0
   real 4     ELON1         ELAT0
   real 5     ELAT1         ELONC
   real 6     ELON2         ELATC
   real 7     ELAT2         EDELX
   real 8     ELON0         EDELY
   real 9     ELAT0         (ELX)
   real 10     ERPK         (ELY)
   real 11   (NSOTRP)      (EXWN)
   real 12   (NGIVO)       (EYWN)
   real 13    (ELX)     ELON1 (West)
   real 14    (ELY)     ELAT1 (South)
   real 15    EDELX     ELON2 (East)
   real 16    EDELY     ELAT2 (North)
   real 17    (EXWN)          –
   real 18    (EYWN)          –

PLONR,PLATR described rotation of the pole. This is obsolete since 2005 or so.

### CADRE\_FOCOHYBRID

This article contains the definition of the hybrid levels used by the
model. Fields on model levels can be recognised by the prefix “S”
followed by the level as a 3 digit integer, for instance
`S010TEMPERATURE` is the temperature at model level 10 (where level 1 is
the top of the atmosphere).

I do not intend to explain the details of hybrid levels in this text.
But they are defined using 2 vectors of length (NLEV+1) and 1 reference
pressure.

  ------------- ----------- --------------------
     real 1      Ref Pres.   Usually 1025000 Pa
   real NLEV+1       A      
   real NLEV+1       B      
  ------------- ----------- --------------------

It is interesting to note that the values of A as coded in the FA frame,
are different from the values that you have to provide in a namelist for
an ALADIN run! They are rescaled by the reference pressure!

### frame name

The name of this record is the “frame name”. This is e.g. used in model
runs to distinguish different post-processing domains. But for offline
decoding of data, it is hardly of any use. The data in this record is 1
(64bit) integer, signifying the version of the FA files. To my
knowledge, this is still 1 in all cases.

### DATE-DES-DONNEES

This is the record that gives date and time information and is almost
identical to the time information in GRIB-1 format. It is a set of 11
integers. The first 5 give the time of the forecast run (dissemination
time).

  -------- -------------------------------------------------------
   int 1   year (including century)
   int 2   month
   int 3   day
   int 4   hour
   int 5   minute
   int 6   forecast range unit (1)
   int 7   P1
   int 8   P2
   int 9   time range indicator
   int 10  Number included in average, when int 9 (Code table 5)
           indicates an average or accumulation;
           otherwise set to zero
   int 11  Number missing from averages or accumulations
  -------- -------------------------------------------------------

The sixth integer gives the unit for the forecast range. Usually, this
is 1 (hours). Other possibilities correspond to those from GRIB code
table 4:

  ----- ----------------------------------
    0   Minute
    1   Hour
    2   Day
    3   Month
    4   Year
    5   Decade (10 years)
    6   Normal (30 years)
    7   Century (100 years)
   8-9  Reserved
   10   3 hours
   11   6 hours
   12   12 hours
   254  Second
        (use DATX-DES-DONNEES in stead!)
  ----- ----------------------------------

The time range indicator is described in GRIB1 code table 5:

  ----- ---------------------------------------------------------------------------------------
    0   fc valid for time + P1 (P1 &gt; 0), or Uninitialized analysis
    1   Initialized analysis product for reference time (P1 $=$ 0)
    2   Product with a valid time ranging between reference time + P1 and reference time + P2
    3   Average (ref. time + P1 to + P2)
    4   Accumulation (ref. time + P1 to + P2) valid at ref. time + P2
    5   Difference (ref. time + P2 minus + P1) valid at ref. time + P2
   10   product valid at ref. time + P1
   ...  ...
  ----- ---------------------------------------------------------------------------------------

There is one difference with code table 5. In GRIB1, the integers are
only 1 octet, so for the case of (10), the forecast range is in fact
P1+256\*P2. This may be useful for longer forecast ranges. However, that
is not necessary here, as the integers are 64 bit.

The last 2 integers are usually 0, except for cumulative fields (fluxes,
precipitation, max/min temperatures over a certain period,...).

int 10 is nowadays set to the time range of the previously created file.
This is important e.g. for fields that give the maximum value since the
previous output.

**NOTE: if the range is not in hours (recent high-resolution may use
shorter output steps) then the hour values are simply rounded and the
correct data is given in `DATX-DES-DONNEES`.**

### DATX-DES-DONNEES

This meta-data field was introduced in **cycle 40**. It gives additional
time information that does not fit in the previous set of 11 integers. I
haven’t found much information yet, but a few fields have been
identified (see `fandai.F90 `):

  -------- ------------------------------------------------------------------------
   int 1   1 (signifies that there is extended data information?)
   int 2   unused (0) ?
   int 3   seconds (additional to Y-M-D H:M already given in `DATE-DES-DONNEES`.)
   int 4   lead time in seconds
   int 5   P1 lead time in seconds (max value since ...) ?
   int 6   P2 lead time in seconds (accumulated since ...) ?
   int 7   timestep in seconds
   int 8   format (0) ?
   int 9   unused (0) ?
   int 10  unused (0) ?
   int 11  unused (0) ?
  -------- ------------------------------------------------------------------------

More information will be added if I get it...

##Data encoding

The data in an FA file may be either in grid point format or spectral.
In either case, the data might be GRIB-compressed or not. For spectral
data there is also the possibility of compressing only the smaller
scales. In all, this gives us a list of different cases to consider.
Some are quite easy, others less so.

The first two entries (64 bit integers) define the kind of data that is
encoded in the field. Integer 1 gives the kind of encoding (none,
internal GRIB or GRIBEX ...). Value 2 is more common than 1. This “FA
extended” GRIB encoding adds separate *min* and *max* values that are
used for more exact reconstruction of the original values. If the GRIB
type is 1, the minimum value is read from the GRIB message itself, and
the maximum is derived from that. The second integer indicates whether
it is spectral or grid point data. Depending on those 6 combinations,
there may be more integers: the number of bits per encoded value
(KNBITS), min and max value in the GRIB encoding, non-compressed scale
(KSTRON).

   entry   gp, no GRIB   gp, GRIB       gp, FA-GRIB       sp, no GRIB   sp, GRIB   sp, FA-GRIB
  ------- ------------- ---------- --------------------- ------------- ---------- -------------
   int 1     $\le 0$        1       $\ge 2$ (usually 2)     $\le 0$        1         $\ge 2$
   int 2        0           0                0                 1           1            1
     3       (data)       KNBITS          KNBITS            (data)       KNBITS      KNBITS
     4                    (data)            min                          KSTRON      KSTRON
     5                                      max                          KPUILB      KPUILB
     6                                    (data)                         (data)        min
     7                                                                                 max
     8                                                                               (data)

**WARNING** in cy43t2, new possibilities have appeared. The first
integer can now have values $-1, 0, 1, 2, 3, 4$. $-1$ and $3$ signify a
different ordering of spectral coefficients. This needs further
investigation. For spectral fields, int2$= -1, 3$ signifies that the
spectral coefficients are ordered differently (i.e. like in the model
vector).

### Uncompressed data

This is by far the simplest case. Especially grid point data is simple:
every value is written as a big endian double. So the data section is
simply a dump of NxM points.

For spectral data, we also have a series of doubles. But now we have to
know exactly how the spectral components are ordered. That is quite
non-trivial. And in fact in rare cases the order may be different (e.g.
in forecast differences as used by FESTAT, the ordering is that of the
model, while in FA files the order is changed).

### GRIB encoding

The GRIB encoding used by the FA library is not exactly the same as in
the WMO documentation. It appears to be GRIB-0 as opposed to GRIB-1 as
used in e.g. GRIBEX. There is an option to use GRIBEX encoding in the FA
file in stead of the old GRIB-0, but I am not aware whether this is ever
used.

In fact only the actual bit-compression of GRIB is used. We are not
interested in the identification section, as the meta-data is in the FA
frame anyway! If the GRIB type is $\ge 2$, even the minimum value and
scaling are not taken from the GRIB part but from the FA frame. The only
data we need is the number of bits per value. That is also given in the
non-compressed header data.

Therefore, we can usually simply skip the **24 octets** at the beginning
and go straight to the actual compressed data section. A little check
(e.g. same number of bits per value?) may be appropriate, but don’t
spend too much time on that. If GRIB type is 1, you do need to read the
minimum value (an IBM-formatted float) and scaling.

**NOTE: in cy43t1, grib-1 and even grib-2 were introduced. I’ll document
that as soon as I found out how it works.**

### Spectral data

If there is no GRIB compression, the spectral components are ordered in
groups of 4 real values that have to be combined to form the complex
Fourier components.

The hardest part is re-ordering the spectral components of a partially
compressed field. First of all: not all components are compressed: only
those with a total wavenumber larger than KSTRON. So first the data is
split into 2 streams.

First comes the compressed part. Before GRIB-encoding, this is
additionally rescaled by a certain power (KPUILB) of the laplacian. Then
it is encoded as a GRIB record. After some padding, the remainder of the
non-compressed fourier components are added.

So to decode this:

1.  At the start of the data section, read the GRIB record and
    decode it.

2.  Rescale the components by their laplacian to some power (i.e.
    multiply by a power of the wavenumbers). For this, you have to be
    careful to have the right N and M value for every
    component (headache).

3.  Merge this data with the uncompressed part (more headache).

4.  Do FFT (as the code takes over, the headache ebbs away).

This has been coded in R.

#ECHKEVO

ECHKEVO files are a bit more obscure. It is possible to write some
fields and grid points to an output file every time step (namelist
`NAMCHK`), with very limited overhead.

This output is gathered in a single file called something like
ICMSHnnnnCHKOUT2.

**NOTE:** at the time of writing, the echkevo implementation is not too
well maintained and, frankly, terribly coded. There is a particularly
nasty bug (feature): meta data is in 32 bit integers, but exported to
file as 64 bit. So two numbers are combined in a single value. That
would not be a problem if all computers were big endian. But a
company named INTEL (you may have heard of it) produces some quite
popular hardware that is little endian. If a model run is performed on
such hardware, the values are written as a byte-swapped 64bit integer,
effectively swapping the values two by two. The worst part of it is that
this is invisible a posteriori. You have to know that if the ECHKEVO
output was *produced* on a little endian machine, the values are swapped
two by two. Nowadays that is in fact the most common case.

Anyway, if you get over the nausea and still want to use ECHKEVO, the
remainder of this section will explain how the data is stored.

The good news is: it’s an LFI-style file.

##Frame

Just like any FA file.

##Info

Detailed information are stored in a data fields called
`REDOCU0000000000`. It’s just a list of (big-endian, possibly swapped)
32bit integers.

`RESTEP0000000000` contains the time step as a single big-endian double.

`REGEOM0000000000` contains some geographical information of the
exported grid points.

##Data

Data fields are called `REnnnn0000000000`, where nnnn is the timestep.
