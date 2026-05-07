#all formulas
red.F <- function(red) { red }
green.F <- function(green) { green }
blue.F <- function(blue) { blue }
RedEdge.F <- function(RedEdge) { RedEdge }
nir.F <- function(nir) { nir }

ndvi.F <- function(nir, red) { (nir - red) / (nir + red) }
# alias because earlier loop used ndvi.F2
#ndvi.F2 <- ndvi.F

ndre.F <- function(nir, RedEdge) { (nir - RedEdge) / (nir + RedEdge) }
dvi.F  <- function(nir, red) { nir - red }
sr.F   <- function(nir, red) { nir / red }

arvi.F <- function(nir, red, blue) { (nir - (2 * red) + blue) / (nir + (2 * red) + blue) }
arvi2.F <- function(nir, red) { -0.18 + 1.17 * ((nir - red) / (nir + red)) }

osavi.F <- function(nir, red) { (1.16 * (nir - red)) / (nir + red + 0.16) }
msavi.F <- function(nir, red) { ((2 * nir + 1) - sqrt((2 * nir + 1)^2 - 8 * (nir - red))) / 2 }
atsavi.F <- function(nir, red) {
  L <- 0.08; a <- 1.22; b <- 0.03
  num <- a * (nir - a * red - b)
  den <- (a * nir + red - a * b + L * (1 + a^2))
  num / den
}

ccci.F <- function(nir, RedEdge, red) { ((nir - RedEdge) / (nir + RedEdge)) / ((nir - red) / (nir + red)) }
cig.F <- function(nir, green) { (nir / green) - 1 }
cire.F <- function(nir, RedEdge) { (nir / RedEdge) - 1 }
ari2.F <- function(nir, green, RedEdge) { (1 / green - 1 / RedEdge) * nir }
cvi.F <- function(nir, red, green) { (nir * red) / (green^2) }
gemi.F <- function(nir, red) {
  part1 <- 2 * (nir^2 - red^2) + 1.5 * nir + 0.5 * red
  frac <- part1 / (nir + red + 0.5)
  gemi_part1 <- frac * (1 - 0.25 * frac)
  gemi_part2 <- (red - 0.125) / (1 - red)
  gemi_part1 - gemi_part2
}

#evi.F <- function(nir, red, blue) { 2.5 * ((nir - red) / (nir + 6 * red - 7.5 * blue + 1)) }
evi2.F <- function(nir, red) { 2.5 * (nir - red) / (nir + 2.4 * red + 1) }
wdrvi.F <- function(nir, red) { (0.1 * nir - red) / (0.1 * nir + red) }

exg.F <- function(green, red, blue) { 2 * green - red - blue }
cive.F <- function(red, green, blue) { 0.441 * red - 0.811 * green + 0.385 * blue + 18.78745 }
gli.F <- function(green, red, blue) { (2 * green - red - blue) / (2 * green + red + blue) }
ri.F <- function(red, green) { (red - green) / (red + green) }
ifindex.F <- function(red, green, blue) { (2 * red - green - blue) / (green - blue) }
intensity.F <- function(red, green, blue) { (1 / 30.5) * (red + green + blue) }

gdvi.F <- function(nir, green) { nir - green }
bndvi.F <- function(nir, blue) { (nir - blue) / (nir + blue) }
gbndvi.F <- function(nir, green, blue) { (nir - (green + blue)) / (nir + (green + blue)) }
grndvi.F <- function(nir, green, red) { (nir - (green + red)) / (nir + (green + red)) }
pndvi.F <- function(nir, green, red, blue) { (nir - (green + red + blue)) / (nir + (green + red + blue)) }
rbndvi.F <- function(nir, red, blue) { (nir - (red + blue)) / (nir + (red + blue)) }
grvi.F <- function(nir, green) { (nir - green) / (nir + green) }
gsavi.F <- function(nir, green) { (1.5 * (nir - green)) / (nir + green + 0.5) }
ndvi_rededge.F <- function(RedEdge, red) { (RedEdge - red) / (RedEdge + red) }

LRVI.F <- function(nir, red) { log10(nir / red) }
rgr.F <- function(red, green) { red / green }
io.F <- function(red, blue) { red / blue }
rri1.F <- function(nir, RedEdge) { nir / RedEdge }
rri2.F <- function(RedEdge, red) { RedEdge / red }
sqrtirr.F <- function(nir, red) { sqrt(nir / red) }
ipvi.F <- function(nir, red) { nir / (nir + red) }

#norms / misc
normG.F <- function(green, nir, red) { green / (nir + red + green) }
normNIR.F <- function(nir, red, green) { nir / (nir + red + green) }
normR.F <- function(red, nir, green) { red / (nir + red + green) }
tndvi.F <- function(nir, red) { sqrt(((nir - red) / (nir + red)) + 0.5) }

# Master wrapper
make_all_vis <- function(red, green, blue, RedEdge, nir) {
  list(
    Red = red.F(red),
    Green = green.F(green),
    Blue = blue.F(blue),
    RedEdge = RedEdge.F(RedEdge),
    NIR = nir.F(nir),
    NDVI = ndvi.F(nir, red),
    NDRE = ndre.F(nir, RedEdge),
    DVI = dvi.F(nir, red),
    SR = sr.F(nir, red),
    ARVI = arvi.F(nir, red, blue),
    ARVI2 = arvi2.F(nir, red),
    OSAVI = osavi.F(nir, red),
    MSAVI = msavi.F(nir, red),
    ATSAVI = atsavi.F(nir, red),
    GSAVI = gsavi.F(nir, green),
    GRVI = grvi.F(nir, green),
    NDVIrededge = ndvi_rededge.F(RedEdge, red),
    CCCI = ccci.F(nir, RedEdge, red),
    CIG = cig.F(nir, green),
    CIRE = cire.F(nir, RedEdge),
    ARI2 = ari2.F(nir, green, RedEdge),
    CVI = cvi.F(nir, red, green),
    GEMI = gemi.F(nir, red),
    #EVI = evi.F(nir, red, blue),
    EVI2 = evi2.F(nir, red),
    WDRVI = wdrvi.F(nir, red),
    ExG = exg.F(green, red, blue),
    CIVE = cive.F(red, green, blue),
    GLI = gli.F(green, red, blue),
    RI = ri.F(red, green),
    IF = ifindex.F(red, green, blue),
    #Hue = hue.F(red, green, blue),
    Intensity = intensity.F(red, green, blue),
    GDVI = gdvi.F(nir, green),
    BNDVI = bndvi.F(nir, blue),
    GBNDVI = gbndvi.F(nir, green, blue),
    GRNDVI = grndvi.F(nir, green, red),
    PNDVI = pndvi.F(nir, green, red, blue),
    RBNDVI = rbndvi.F(nir, red, blue),
    LRVI = LRVI.F(nir, red),
    RGR = rgr.F(red, green),
    IO = io.F(red, blue),
    RRI1 = rri1.F(nir, RedEdge),
    RRI2 = rri2.F(RedEdge, red),
    SQRTIRR = sqrtirr.F(nir, red),
    IPVI = ipvi.F(nir, red),
    #GARI = gari.F(nir, green, blue, red),
    NormG = normG.F(green, nir, red),
    NormNIR = normNIR.F(nir, red, green),
    NormR = normR.F(red, nir, green),
    #CI = ci.F(red, blue),
    TNDVI = tndvi.F(nir, red)
  )
}