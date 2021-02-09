require(ggplot2)
require(ggthemes)
require(ggstance)
require(extrafont)
require(complexplus)
require(OpenImageR)
require(envi)
require(float)
library(caTools)
library(raster)
library(rgdal)
library(rasterVis)
library(RColorBrewer)

#font_import()
#loadfonts()
##---------------------------------------------------------------------------------------
readComplexFloat <- function(f, size, endian){
  X2 <- readBin(f, double(), n = 2*size, size = 4, endian = endian) # caso 4
  X = complex(length.out = size)
  
  for(i in 1:size) {
    twoi = 2*i
    X[i] = complex(real = X2[twoi-1], imaginary=X2[twoi])
  }
  
  return(X)
}

myread.ENVI <- function (filename, headerfile = paste(filename, ".hdr", sep = ""))
{
  nCol <- nRow <- nBand <- data.type <- header.offset <- byte.order <- (-1)
  interleave = "bsq"
  if (!file.exists(headerfile))
    stop("read.ENVI: Could not open input header file: ",
         headerfile)
  Lines = read.table(headerfile, sep = "=", strip.white = TRUE,
                     row.names = NULL, as.is = TRUE, fill = TRUE)
  Fields = c("samples", "lines", "bands", "data type", "header offset",
             "interleave", "byte order")
  for (i in 1:nrow(Lines)) {
    Lab = tolower(Lines[i, 1])
    Lab = gsub("[ ]+", " ", Lab)
    j = match(Lab, Fields)
    Val = Lines[i, 2]
    if (length(j) == 1)
      switch(j, nCol <- as.integer(Val), nRow <- as.integer(Val),
             nBand <- as.integer(Val), data.type <- as.integer(Val),
             header.offset <- as.integer(Val), interleave <- gsub(" ",
                                                                  "", Val), byte.order <- as.integer(Val))
  }
  if (nCol <= 0 | nRow <= 0 | nBand <= 0)
    stop("read.ENVI: Error in input header file ", headerfile,
         " data sizes missing or incorrect", nRow, nCol, nBand)
  if (!(data.type %in% c(1, 2, 3, 4, 5, 6, 9, 12)))
    stop("read.ENVI: Error in input header file ", headerfile,
         " data type is missing, incorrect or unsupported ")
  ieee = if (.Platform$endian == "big")
    1
  else 0
  endian = if (ieee == byte.order | byte.order < 0)
    .Platform$endian
  else "swap"
  size = nRow * nCol * nBand
  if (!file.exists(filename))
    stop("read.ENVI: Could not open input file: ", filename)
  f = file(filename, "rb")
  if (header.offset > 0)
    readBin(f, raw(), n = header.offset)
  switch(data.type,
         X <- readBin(f, integer(), n = size, size = 1, signed = FALSE), # caso 1
         X <- readBin(f, integer(), n = size, size = 2, endian = endian), # caso 2
         X <- readBin(f, integer(), n = size, endian = endian), # caso 3
         X <- readBin(f, double(), n = size, size = 4, endian = endian), # caso 4
         X <- readBin(f, double(), n = size, endian = endian), # caso 5
         X <- readComplexFloat(f, size = size, endian = endian), # caso 6
         ,  # caso 7
         ,  # caso 8
         X <- readBin(f, complex(), n = size, endian = endian), # caso 9
         , # caso 10
         , # caso 11
         X <- readBin(f, integer(), n = size, size = 2, endian = endian, signed = FALSE) # caso 12
  )
  close(f)
  Fields = c("bil", "bip", "bsq")
  j = match(interleave, Fields)
  if (length(j) == 0)
    stop("read.ENVI: Error in input header file ", headerfile,
         " incorrect interleave type")
  switch(j, {
    dim(X) <- c(nCol, nBand, nRow)
    X <- aperm(X, c(3, 1, 2))
  }, {
    dim(X) <- c(nBand, nCol, nRow)
    X <- aperm(X, c(3, 2, 1))
  }, {
    dim(X) <- c(nCol, nRow, nBand)
    X <- aperm(X, c(2, 1, 3))
  })
  if (nBand == 1)
    dim(X) = c(nRow, nCol)
  return(X)
}

para01 <- function(x) {
  valores <- range(x)
  y = (x - valores[1]) / (valores[2] - valores[1])
  y
}


##----------------------------------------------------------------------------
main_dir = 'D:\\MF4CF\\T3\\'

# Read coherence elements (T3)
T11 <- myread.ENVI(paste0(main_dir,"T11.bin"), 
                   paste0(main_dir,"T11.bin.hdr"))
T12_imag <- myread.ENVI(paste0(main_dir,"T12_imag.bin"), 
                        paste0(main_dir,"T12_imag.bin.hdr"))
T12_real <- myread.ENVI(paste0(main_dir,"T12_real.bin"), 
                        paste0(main_dir,"T12_real.bin.hdr"))

# T12 = T12_real + T12_imag*1i
# T21 = Conj(T12)

T13_imag <- myread.ENVI(paste0(main_dir,"T13_imag.bin"), 
                        paste0(main_dir,"T13_imag.bin.hdr"))
T13_real <- myread.ENVI(paste0(main_dir,"T13_real.bin"), 
                        paste0(main_dir,"T13_real.bin.hdr"))

# T13 = T13_real + T13_imag*1i
# T31 = Conj(T13)

T22 <- myread.ENVI(paste0(main_dir,"T22.bin"), 
                   paste0(main_dir,"T22.bin.hdr"))
T23_imag <- myread.ENVI(paste0(main_dir,"T23_imag.bin"), 
                        paste0(main_dir,"T23_imag.bin.hdr"))
T23_real <- myread.ENVI(paste0(main_dir,"T23_real.bin"), 
                        paste0(main_dir,"T23_real.bin.hdr"))

# T23 = T23_real + T23_imag*1i
# T32 = Conj(T23)

T33 <- myread.ENVI(paste0(main_dir,"T33.bin"), 
                   paste0(main_dir,"T33.bin.hdr"))

Nrw = nrow(T11)
Ncl = ncol(T11)

##----------------------------------------------------------------------------
# selection of window
wsi = 7; # window size

kernel = matrix(1, nrow = wsi, ncol = wsi) / (wsi*wsi)


t11s = convolution(T11, kernel, "same")
T12_imag_filt = convolution(T12_imag, kernel, "same")
T12_real_filt = convolution(T12_real, kernel, "same")

t12s = T12_real_filt + T12_imag_filt*1i
t21s = Conj(t12s)

T13_imag_filt = convolution(T13_imag, kernel, "same")
T13_real_filt = convolution(T13_real, kernel, "same")

t13s = T13_real_filt + T13_imag_filt*1i
t31s = Conj(t13s)

t22s = convolution(T22, kernel, "same")

T23_imag_filt = convolution(T23_imag, kernel, "same")
T23_real_filt = convolution(T23_real, kernel, "same")

t23s = T23_real_filt + T23_imag_filt*1i
t32s = Conj(t23s)

t33s = convolution(T33, kernel, "same")

det_T3 = t11s*(t22s*t33s-t23s*t32s)-t12s*(t21s*t33s-t23s*t31s)+t13s*(t21s*t32s-t22s*t31s)
trace_T3 = t11s + t22s + t33s
m1 = Re(sqrt(1-(27*(det_T3/(trace_T3^3)))))

k11_f = (t11s + t22s + t33s)/2
k44_f = (-t11s + t22s + t33s)/2
k14_f = Im(t23s)

s0_f = trace_T3
dop_f = m1

val1 = (4*dop_f*k11_f*k44_f)/(k44_f^2 - (1 + 4*dop_f^2)*k11_f^2)
val2 = abs(k14_f)/(k11_f)

theta_f = atan(val1)*(180/pi) # separation for surface and dbl
tau_f = atan(val2)*(180/pi) # separation for helix

pc_f = dop_f*s0_f*(sin(2*tau_f*(pi/180)))
pv_f = (1-dop_f)*s0_f
res_pow = s0_f - (pc_f + pv_f)
ps_f = (res_pow/2)*(1+sin((2*theta_f*(pi/180))))
pd_f = (res_pow/2)*(1-sin((2*theta_f*(pi/180))))

#-----------------------------------------------------------------------------
# Save in envi format
dir.create(file.path(main_dir, 'MF4CF_Rext'), showWarnings = FALSE)
save_dir = paste0(main_dir,'MF4CF_Rext\\')

r <- raster(theta_f)
r <- writeRaster(r, filename=paste0(save_dir,'theta_f.envi'), overwrite=TRUE)

r <- raster(tau_f)
r <- writeRaster(r, filename=paste0(save_dir,'tau_f.envi'), overwrite=TRUE)

r <- raster(pc_f)
r <- writeRaster(r, filename=paste0(save_dir,'pc_f.envi'), overwrite=TRUE)

r <- raster(pv_f)
r <- writeRaster(r, filename=paste0(save_dir,'pv_f.envi'), overwrite=TRUE)

r <- raster(ps_f)
r <- writeRaster(r, filename=paste0(save_dir,'ps_f.envi'), overwrite=TRUE)

r <- raster(pd_f)
r <- writeRaster(r, filename=paste0(save_dir,'pd_f.envi'), overwrite=TRUE)