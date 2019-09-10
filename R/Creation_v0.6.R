#' Package for creating simulated ascii land covers
#'
#' This package allows for presicie control of generation of focal habitats
#' as well as filling in the matrix habitats inbetween the focal patches.
#'
#'
#' @details
#'     The main function i is \link[PhiCor]{CorIndex.all.plusP}.
#'     This will calculate all phi corelation of association values and
#'     purmutate to give p-values.
#'
#'     The other Phi coefficent functions are useful in testing or calculating only part
#'     of the full outputs from \link[PhiCor]{CorIndex.all.plusP}:
#'     \itemize{
#'         \item \link[PhiCor]{CorIndex} calculates the inputs that don't change between habitats.
#'         \item \link[PhiCor]{CorIndex.TargetVar} calculates the additional input specific to a habitat.
#'         \item \link[PhiCor]{CorIndex.groupEqual} calculates the group-equalised phi for a habitat.
#'         \item \link[PhiCor]{CorIndex.nonEqual} calculates the non group-equalised phi for a habitat.
#'         \item \link[PhiCor]{CorIndex.xPerm} returns the null phi distribution for a habitat.
#'         \item \link[PhiCor]{CorIndex.all} returns the non and group-equalised phi values for all habitats
#'             without p-values and is therefore quicker.
#'     }
#'
#'
#'     The functions \link[PhiCor]{lookup} and \link[PhiCor]{FilePathInTree}
#'     might be used elsewhere.
#' @author Jordan Chetcuti
#' @references
#'
#' @references
#'
#' @keywords internal
#'
"_PACKAGE"
#> [1] "_PACKAGE"


#' header is the header for a asc file. This is a bit of a bodge, but can be added to
#' an array to give and ascii file.
#'
#'
#' @format Six character objects defining an ascii
#' @usage The data is automatically used and has no real independent use.
"header"



#' ArrayToAsc converts an array to asc format modifying and stitching a header on.
#'
#'This adds on a simple header to an array to create an asc spatial dataset and saves
#'this to the hard drive.
#'
#' @param inArray array - an input array that will provide cell values for the asc file.
#' @param filename character - a string giving the filename and path.
#' @param overwrite boolean - by default set as NULL. If this is set to true, then the file
#' can overwrite another ascii file with same name.
#' @param runno boolean - as part of the process a tempory file is created. If multiple
#' ascii files are being created simultanously then an identifier may stop a conflict of temp
#' files.
#' @return character - currently returns the filepath
#' @examples
#'
#' TBA...
#'
#' @export
ArrayToAsc <- function(inArray, filename, overwrite=FALSE, runno = NULL){

  if(!is.null(runno)){
    TempFielName = paste("test",runno, ".asc", sep="")
  }else{
    TempFielName = "test.asc"
  }


  write.table(inArray, TempFielName, sep = " ", row.names = FALSE, col.names = FALSE)
  body = readLines(TempFielName)
  file.remove(TempFielName)
  header = dget("header")
  header <- gsub("NROWS 1020", paste("NROWS", nrow(inArray), sep = " "))
  header <- gsub("NCOLS 1020", paste("NCOLS", ncol(inArray), sep = " "))
  OutPut = c(header,body)

  if(overwrite == TRUE & file.exists(filename)){
    #Delete file if it exists
    file.remove(filename)
  }


  writeLines(OutPut, filename)
  return(filename)
}


#' ExtractWindow extracts a defined region of an array.
#'
#' This extracts a number of cells from an array. Typically used in a moving window for example.
#'
#' @param InMatrix array - an input array from which a subset will be taken.
#' @param y integer - central y location for subset.
#' @param x integer - central x location for subset.
#' @param ySize integer - height of window
#' @param xSize integer - width of window
#' @return array - a array smaller than the input array of the 3 by 3 by default or larger if specified.
#' @examples
#'
#' TBA...
#'
#' @export
ExtractWindow = function(InMatrix, y, x, ySize = 3, xSize = 3){
  xExtent = ncol(InMatrix)
  yExtent = nrow(InMatrix)


  OutMatrix = matrix(0,ySize,xSize, byrow=TRUE)

  xSize = (xSize - 1) / 2
  ySize = (ySize - 1) / 2
  countX = 1
  for(xRange in (x-xSize):(x+xSize)){
    countY = 1
    for(yRange in (y-ySize):(y+ySize)){
      if(xRange>0 && xRange <= xExtent && yRange>0 && yRange <= yExtent){
        OutMatrix[ countY, countX] = InMatrix[yRange, xRange]
      }
      countY = countY + 1
    }
    countX = countX + 1
  }
  return(OutMatrix)
}

#' GrowBuffer adds a cell of habitat when growing habitat patches.
#'
#' I wrote this and I've lost track of exactly what it does...
#'
#' @param yLocal array - an input array from which a subset will be taken.
#' @param xLocal integer - central y location for subset.
#' @param seed integer - central x location for subset.
#' @param Extent integer - height of window
#' @param World integer - width of window
#' @param GrowWorld array
#' @param WeightWorld array
#' @return list of arrays - the list contains two arrays, the updated GrowWorld and WeightWorld.
#' @examples
#'
#' TBA...
#'
#' @export
GrowBuffer <- function(yLocal, xLocal, seed, Extent, World, GrowWorld, WeightWorld){

  if(xLocal==1){
    minx = xLocal
  }else if(xLocal>1){
    minx = xLocal-1
  }

  if(xLocal==Extent){
    maxx = xLocal
  }else if(xLocal<Extent){
    maxx = xLocal+1
  }

  if(yLocal==1){
    miny = yLocal
  }else if(yLocal>1){
    miny = yLocal-1
  }

  if(yLocal==Extent){
    maxy = yLocal
  }else if(yLocal<Extent){
    maxy = yLocal+1
  }





  xRange = minx:maxx

  #cdef vector[int] yRange
  yRange = miny:maxy

  for (xValue in xRange){
    for (yValue in yRange){
      if (xValue>0 && xValue<=Extent && yValue>0 && yValue<=Extent){
        if(World[yValue, xValue]==0 && (!(xValue == xLocal && yValue == yLocal)) && GrowWorld[yValue,xValue]==0){
          GrowWorld[yValue,xValue] <- seed
          WeightWorld[yValue,xValue] <- WeightWorld[yValue,xValue] + 1
        }else if (World[yValue,xValue]==0 && (!(xValue == xLocal && yValue == yLocal))){
          WeightWorld[yValue,xValue] <- WeightWorld[yValue,xValue] + 1
        }
      }
    }
  }

  #print(plot(raster(GrowWorld)))
  outList <- list("GrowWorld" = GrowWorld, "WeightWorld" = WeightWorld )

  return(outList)
}

AddAllSeedPoints <- function(AllyLocal, AllxLocal, World, GrowWorld, WeightWorld, Extent){
  #seed, nonseed, Locals, World, GrowDict, WeightWorld, Extent):
  #print(length(AllyLocal))
  RangeSeeds = 1:length(AllyLocal)
  for(seed in RangeSeeds){
    #print(seed)
    nonseed = RangeSeeds[RangeSeeds !=seed ]
    yLocal = AllyLocal[seed]
    xLocal = AllxLocal[seed]
    while(World[ yLocal, xLocal]>0 | GrowWorld[ yLocal, xLocal]>0 | sum(ExtractWindow(World, yLocal, xLocal) %in% nonseed)>0){
      if(xLocal==1){
        xLocal <- xLocal + 1
      }else if(xLocal==Extent){
        xLocal <- xLocal - 1
      }else{
        xLocal <- sample(c(-1,0, 1), 1) + xLocal
      }

      if(yLocal==1){
        yLocal <- yLocal + 1
      }else if(yLocal==Extent){
        yLocal <- yLocal - 1
      }else{
        yLocal <- sample(c(-1,0,1), 1) + yLocal
      }
    }
    World[yLocal, xLocal] <- seed
    Recieved <- GrowBuffer(yLocal, xLocal, seed, Extent, World, GrowWorld, WeightWorld)
    GrowWorld <- Recieved$GrowWorld
    WeightWorld <- Recieved$WeightWorld
  }
  #plot(raster(GrowWorld))
  Recieved$World <- World
  return(Recieved)
}

#' @import mobsim
#' @import spatstat
#' @export
CreatSeeds <- function(PatchDistribution, SeedPoints,Extent){
  # library(mobsim)
  # library(spatstat)
  if(is.null(PatchDistribution)||PatchDistribution=="Random"||PatchDistribution=="Poisson"){
    #print("Random")
    SeedXY = sim_poisson_coords(rep(1,SeedPoints), xrange = c(1,Extent), yrange = c(1,Extent))$`census`
  }else if(PatchDistribution=="Clustered"|PatchDistribution=="Thomas"){
    #print("Clustered")
    "Currently I have this set up to make the extremes of the number of clusters
    less likely. We might need to change this."
    Mother = sample(1:SeedPoints, 1, prob = dnorm(1:SeedPoints, mean=(0.5 + (SeedPoints/2)), sd=(SeedPoints/2)))
    #print(Mother)
    SeedXY <- sim_thomas_coords(c(SeedPoints), mother_points = Mother, sigma = 0.02,xrange = c(1,Extent), yrange = c(1,Extent))$`census`
  }else if(PatchDistribution=="Regular"|PatchDistribution=="Strauss"|PatchDistribution=="Dispersed"){
    "I'm not sure what effect the beta value really has, therefore I have left it as 200.
    I'm using 1/sqrt(number of seed points) from May et al. (2019). I've used a beta distribution
    for the gamma value. This value if zero is totally distributed, and if it is one then the
    points are randomly placed again. This beta distribution should mean that the value is likely
    to be closer to zero than one most of the time."
    nr = 10000 # no. of repetitions in MCMC run
    nv = 0
    #p2 = list(cif="strauss",par=list(beta=200,gamma=rbeta(1, 1, 3),r=1/(SeedPoints^0.5)),w=c(1,Extent,1,Extent))
    "Changed my mind, I will force them to be regular."
    #I calculated an adjustement to the way of predicting R from May et al.(2019)
    #This adjustment maximises the Nearest Neighbour Index so that
    RVal = 0.769019 * (Extent/(SeedPoints^0.5))- 0.721628
    p2 = list(cif="strauss",par=list(beta=200,gamma=0.1,r=RVal),w=c(1,Extent,1,Extent))
    sink("file")
    SeedXY <- unclass(rmh(model=p2,start=list(n.start=SeedPoints),control=list(p=1,nrep=nr,nverb=nv)))
    sink()
  }
  return(SeedXY)
}



#' CreateWorld sets up an initial matrix with the focal patches in
#'
#' creates a landscape with zeros in and then adds patches
#'
#' @param SeedPoints integer - number of patches
#' @param PercentageCover float between 0.0 & 1.0 giving the proportion of array that should be patches.
#' @param Extent integer - size of square array, defaults to 250 (large arrays take longer to create.)
#' @param SeedDistribution character - how seed points are distributed to grow the patches. Options are "Random" or "Poisson", "Clustered" or "Thomas", or "Regular", "Strauss" or"Dispersed"
#' and defaults to random.
#' @param PatchDistribution  character -a distribution from which to pick the size of patches from. Opetions are
#' "Gaussian" or"Normal", "Pareto", or the default "Uniform". Uniform gives patches of almost any size and the maximum variation.
#' @param SecondPercentageCover float between 0.0 & 1.0 giving the proportion of array that should be patches.If > 0 creates an additional array based on the seeds of the first with a different proportion.
#' @param normalSD default = 0.01,
#' @param paretoVal defulat = 1
#' @return list of arrays - the list contains two arrays, the updated GrowWorld and WeightWorld.
#' @examples
#'
#' TBA...
#' @import truncnorm
#'
#' @export
CreateWorld <- function(SeedPoints, PercentageCover, Extent, SeedDistribution = NULL, PatchDistribution=NULL, SecondPercentageCover = 0.0, normalSD=0.01, paretoVal = 1){
  # library(truncnorm)

  NumberCells = Extent^2
  RangeSeeds = 1:SeedPoints
  "This allows a distribution to be added for how the patches are grown. We can add different
  distributions here if we need them."
  if(is.null(SeedDistribution)){
    SeedDist = rep(1, SeedPoints)
  }else if(SeedDistribution =="Gaussian" | SeedDistribution =="Normal"){
    #print("NORMAL!!!")
    #SeedDist = rnorm(SeedPoints, mean =  SeedPoints/2, sd=normalSD)
    #SeedDist = SeedDist/sum(SeedDist)

    SeedDist = 100*rtruncnorm(1:SeedPoints, a=0, b=2, mean = 1, sd = normalSD)
    #SeedDist = rnorm(1:SeedPoints, mean=1, sd=normalSD)
  }else if(SeedDistribution =="Pareto"){
    SeedDist = rpareto(1:SeedPoints, paretoVal, 100)
  }else if(SeedDistribution =="Uniform"){
    SeedDist = runif(1:SeedPoints, min=0, max=100)
  }
  #print(SeedDist)
  #The landcover can be created for a pair of land covers with different cover amounts.
  #This if statement means they don't have to be in order of larger and smaller.
  if(PercentageCover>SecondPercentageCover && !PercentageCover==SecondPercentageCover){
    CoverCells = NumberCells * PercentageCover
    SmallerCoverCells = NumberCells * SecondPercentageCover
  }else if(PercentageCover<SecondPercentageCover && !PercentageCover==SecondPercentageCover){
    CoverCells = NumberCells * SecondPercentageCover
    SmallerCoverCells = NumberCells * PercentageCover
  }
  pb <- txtProgressBar(min = 0, max = CoverCells, style = 3)
  #Initiallising the matricies.
  World = matrix(0, nrow = Extent, ncol = Extent,  byrow=TRUE)
  GrowWorld = matrix(0, nrow = Extent, ncol = Extent, byrow=TRUE)
  WeightWorld = matrix(0, nrow = Extent, ncol = Extent, byrow=TRUE)
  Locals = 1:Extent

  SeedXY = CreatSeeds(PatchDistribution, SeedPoints,Extent)

  #print(as.integer(round(SeedXY$y)))
  Recieved <- AddAllSeedPoints(as.integer(round(SeedXY$y)), as.integer(round(SeedXY$x)), World, GrowWorld, WeightWorld, Extent)
  World <- Recieved$World
  GrowWorld <- Recieved$GrowWorld
  WeightWorld <- Recieved$WeightWorld
  #print(plot(raster(GrowWorld)))
  #This grows the patches.
  while(length(World[World > 0])<CoverCells){
    "Randomly choose a patch to grow. If we switch the distribution, then we would get
    some grown more often than others."
    #print(length(World[World > 0]))
    setTxtProgressBar(pb, length(World[World > 0]))
    #This should pick patches that can be grown.
    Possiblities = c()
    while(length(Possiblities)==0){
      targetPatch = sample(RangeSeeds, 1, prob = SeedDist)
      Possiblities = which(GrowWorld==targetPatch , arr.ind=TRUE)
    }
    other = RangeSeeds[RangeSeeds !=targetPatch ]
    "If we change the distribution here, it might grow the patch in less uniform ways.
    This needs some experimentation though. Currently this also adds a weight to fill
    in pixels within a patch."
    ##############################
    #NumberPossibilities = nrow(Possiblities)
    #ATest = rgamma(NumberPossibilities, shape = 1)
    ##############################
    #This did nothing different when multiplying the existing probability and adding it to the exception.
    #1 over the probability that currently fills in the patches, created "frothy" patches. We could have
    #Filled in, not filled in and frothy patches, but I can't see what we would want this. I think leave
    #this as it is.
    APossibileRow = try(sample(1:nrow(Possiblities),1, prob = (WeightWorld[Possiblities]/max(WeightWorld))))
    if(class(APossibileRow)== "try-error"){
      APossibileRow = sample(1:nrow(Possiblities),1)
    }
    APossibility = Possiblities[APossibileRow,]
    if(World[APossibility[1], APossibility[2]]==0 && sum(ExtractWindow(World, APossibility[1], APossibility[2]) %in% other)==0){
      World[APossibility[1], APossibility[2]]=targetPatch
      #yLocal, xLocal, seed, Extent, World, GrowWorld, WeightWorld
      Recieved <- GrowBuffer(APossibility[1], APossibility[2], targetPatch, Extent, World, GrowWorld, WeightWorld)
      GrowWorld <- Recieved$GrowWorld
      WeightWorld <- Recieved$WeightWorld
    }
    GrowWorld[APossibility[1], APossibility[2]]= 0
    if (length(World[World > 0]) == SmallerCoverCells){
      SmallerWorld = World
    }
  }
  if(SmallerCoverCells>0){
    OutList = list("LargerWorld"= World, "SmallerWorld" = SmallerWorld, "GrowWorld" = GrowWorld)
  }else{
    OutList = list("LargerWorld"= World, "SmallerWorld" = NULL, "GrowWorld" = GrowWorld)
  }
  return(OutList)
}

GetDistance <-function(InSeedxy,AvaiableCells){
  #print(InSeedxy)
  if(is.atomic(InSeedxy)){
    SeedX = as.numeric(InSeedxy[1])
    SeedY = as.numeric(InSeedxy[2])
  }else{
    SeedX = as.numeric(InSeedxy$x)
    SeedY = as.numeric(InSeedxy$y)
  }
  #print(SeedX)
  #print(SeedY)
  distance=((SeedX - AvaiableCells[,2])^2 + (SeedY - AvaiableCells[,1])^2)^0.5
  Row = as.numeric(AvaiableCells[which(distance==min(distance), arr.ind=TRUE),])
  return(Row)
}

AddMoreSeedPoints <- function(AllyLocal, AllxLocal, World, GrowWorld, WeightWorld, Extent){
  #seed, nonseed, Locals, World, GrowDict, WeightWorld, Extent):
  #print(length(AllyLocal))
  StartMax = max(World)

  RangeOfPoints = 1:length(AllyLocal)
  RangeSeeds = RangeOfPoints + StartMax
  for(index in RangeOfPoints){
    #print(seed)
    seed = index + StartMax
    nonseed = RangeSeeds[RangeSeeds !=seed ]
    yLocal = AllyLocal[index]
    xLocal = AllxLocal[index]
    while(World[ yLocal, xLocal]>0 | GrowWorld[ yLocal, xLocal]>0 | sum(ExtractWindow(World, yLocal, xLocal) %in% nonseed)>0){
      if(xLocal==1){
        xLocal <- xLocal + 1
      }else if(xLocal==Extent){
        xLocal <- xLocal - 1
      }else{
        xLocal <- sample(c(-1,0, 1), 1) + xLocal
      }

      if(yLocal==1){
        yLocal <- yLocal + 1
      }else if(yLocal==Extent){
        yLocal <- yLocal - 1
      }else{
        yLocal <- sample(c(-1,0,1), 1) + yLocal
      }
    }
    World[yLocal, xLocal] <- seed
    Recieved <- GrowBuffer(yLocal, xLocal, seed, Extent, World, GrowWorld, WeightWorld)
    GrowWorld <- Recieved$GrowWorld
    WeightWorld <- Recieved$WeightWorld
  }
  #plot(raster(GrowWorld))
  Recieved$World <- World
  return(Recieved)
}

#' @import truncnorm
GrowAllBuffers <-function(World,GrowWorld, WeightWorld, CoverCells,RangeSeeds , SeedDistribution=NULL, normalSD=0.1){
  # library(truncnorm)
  Extent<- nrow(World)
  SeedPoints = length(RangeSeeds)
  "This allows a distribution to be added for how the patches are grown. We can add different
  distributions here if we need them."
  if(is.null(SeedDistribution)){
    SeedDist = rep(1, SeedPoints)
  }else if(SeedDistribution =="Gaussian" | SeedDistribution =="Normal"){
    #print("NORMAL!!!")
    SeedDist = rnorm(SeedPoints, mean =  SeedPoints/2, sd=normalSD)
    SeedDist = SeedDist/sum(SeedDist)
  }else if(SeedDistribution =="Uniform"){
    SeedDist = runif(1:SeedPoints, min=0, max=100)
  }

  while(length(World[World > 0])<CoverCells & length(GrowWorld[GrowWorld>0])>0){
    "Randomly choose a patch to grow. If we switch the distribution, then we would get
    some grown more often than others."
    #print(length(World[World > 0]))

    #This should pick patches that can be grown.
    Possiblities = c()
    while(length(Possiblities)==0){
      if(length(RangeSeeds)==length(SeedDist)){
        targetPatch = sample(RangeSeeds, 1, prob = SeedDist)
      }else{
        save(RangeSeeds, paste("RS",format(Sys.time(), "%H:%M:%S"),".rdata", sep = ""))
        save(SeedDist, paste("SD",format(Sys.time(), "%H:%M:%S"),".rdata", sep = ""))
        #targetPatch = sample(RangeSeeds, 1)
      }

      Possiblities = which(GrowWorld==targetPatch , arr.ind=TRUE)
    }
    other = RangeSeeds[RangeSeeds !=targetPatch ]
    "If we change the distribution here, it might grow the patch in less uniform ways.
    This needs some experimentation though. Currently this also adds a weight to fill
    in pixels within a patch."
    ##############################
    #NumberPossibilities = nrow(Possiblities)
    #ATest = rgamma(NumberPossibilities, shape = 1)
    ##############################
    #This did nothing different when multiplying the existing probability and adding it to the exception.
    #1 over the probability that currently fills in the patches, created "frothy" patches. We could have
    #Filled in, not filled in and frothy patches, but I can't see what we would want this. I think leave
    #this as it is.
    APossibileRow = try(sample(1:nrow(Possiblities),1, prob = (WeightWorld[Possiblities]/max(WeightWorld))))
    if(class(APossibileRow)== "try-error"){
      APossibileRow = sample(1:nrow(Possiblities),1)
    }
    APossibility = Possiblities[APossibileRow,]
    if(World[APossibility[1], APossibility[2]]==0 && sum(ExtractWindow(World, APossibility[1], APossibility[2]) %in% other)==0){
      World[APossibility[1], APossibility[2]]=targetPatch
      #yLocal, xLocal, seed, Extent, World, GrowWorld, WeightWorld
      Recieved <- GrowBuffer(APossibility[1], APossibility[2], targetPatch, Extent, World, GrowWorld, WeightWorld)
      GrowWorld <- Recieved$GrowWorld
      WeightWorld <- Recieved$WeightWorld
    }
    GrowWorld[APossibility[1], APossibility[2]]= 0
  }
  OutList = list("World"= World, "SmallerWorld" = NULL, "GrowWorld" = GrowWorld)

  return(OutList)
}

CreateAdditionalWorld <- function(World, LCV, CoverOfPatch, SeedPoints, LCVHabitat, PatchDistribution=NULL,SeedDistribution = "Uniform"){

  "Having created the initial habitat patches of interest these need to be assigned
  a landcover (In a new matrix) and then we need to grow the other habitats.
  I'm going to do this for a single set of matrix variables initially, then we can
  duplicate."
  #This would need to be done for the smaller not the larger if there were two.

  #print(LCV)

  UsedCells <- which(World>0, arr.ind=TRUE)
  MatrixCells <- which(World==0, arr.ind=TRUE)
  "Randomly seeding shouldn't be difficult, just need to randomly choose from this list.
  The other distributions might be a bit trickier. Unless I use the true false matrix
  and then just move a seed point to the nearest true location."
  #Five other habs to start with.

  #How much matrix do we have to work with?
  UsedSize = nrow(UsedCells)


  Extent =nrow(World)

  SeedXY <- CreatSeeds(PatchDistribution, SeedPoints,Extent)
  #if I could find the nearest MatrixCells points to these, then I could move them to these.
  #print(SeedXY)
  SeedXY = as.data.frame(SeedXY[c("x","y")])
  #NewXY = apply(SeedXY,1, f)


  NewXY = apply(SeedXY,1, GetDistance, AvaiableCells=MatrixCells)

  rm(MatrixCells)
  SeedXY$NewY = NewXY[1,]
  SeedXY$NewX = NewXY[2,]
  GrowWorld = matrix(0, nrow = Extent, ncol = Extent, byrow=TRUE)
  WeightWorld = matrix(0, nrow = Extent, ncol = Extent, byrow=TRUE)

  ExistingPatches = max(World)

  Recieved <- AddMoreSeedPoints(SeedXY$NewY,SeedXY$NewX, World, GrowWorld, WeightWorld, Extent)
  World <- Recieved$World
  GrowWorld <- Recieved$GrowWorld
  WeightWorld <- Recieved$WeightWorld
  TotCover = CoverOfPatch + UsedSize
  #print(plot(raster(GrowWorld)))
  RangeNewPatches = 1:SeedPoints + ExistingPatches

  InterestPatches = GrowAllBuffers(World, GrowWorld, WeightWorld, TotCover, RangeNewPatches, SeedDistribution)

  World = InterestPatches$World


  #This adds the second landcover
  LCV = ifelse(World>0 & LCV==0,LCVHabitat ,LCV)

  #print(plot(raster(InterestLCV)))


  OutList = list("Patches"= World, "LCV" = LCV)
  return(OutList)

}

FillBuffer <- function(yLocal, xLocal, seed, Extent, World, GrowWorld){
  if(xLocal==1){
    minx = xLocal
  }else if(xLocal>1){
    minx = xLocal-1
  }

  if(xLocal==Extent){
    maxx = xLocal
  }else if(xLocal<Extent){
    maxx = xLocal+1
  }

  if(yLocal==1){
    miny = yLocal
  }else if(yLocal>1){
    miny = yLocal-1
  }

  if(yLocal==Extent){
    maxy = yLocal
  }else if(yLocal<Extent){
    maxy = yLocal+1
  }
  xRange = minx:maxx
  yRange = miny:maxy
  for (xValue in xRange){
    for (yValue in yRange){
      if (xValue>0 && xValue<=Extent && yValue>0 && yValue<=Extent){
        if(World[yValue, xValue]==0 && (!(xValue == xLocal && yValue == yLocal)) && GrowWorld[yValue,xValue]==0){
          GrowWorld[yValue,xValue] <- seed
        }
      }
    }
  }
  return(GrowWorld)
}

FillZeroPatch <-function(World,GrowWorld, WeightWorld){
  Extent<- nrow(World)

  targetPatch = max(World)

  while(length(GrowWorld[GrowWorld > 0])>0){
    #This should pick patches that can be grown.
    Possiblities = which(GrowWorld>0, arr.ind=TRUE)
    APossibileRow = sample(1:nrow(Possiblities),1)
    APossibility = Possiblities[APossibileRow,]
    GrowWorld  <- FillBuffer(APossibility[1], APossibility[2], targetPatch, Extent, World, GrowWorld)
    World[APossibility[1], APossibility[2]]=targetPatch
    GrowWorld[APossibility[1], APossibility[2]]= 0
  }
  OutList = list("Patches"= World)
  return(OutList)
}

FillAllZeroPatches <-function(Patches){
  Extent<- nrow(Patches)
  Possiblities = which(Patches==0 , arr.ind=TRUE)


  while(nrow(Possiblities)>0){

    CurrentMaxPatch = max(Patches)

    RowChosen = sample(1:nrow(Possiblities),1)

    StartingSeed = Possiblities[RowChosen,]
    #This is basically a seed point. I should then be able to use the same growing code to fill it in.
    #Each one will then be reassesed.
    #Except it will need to be until length(GrowWorld) ==0 (or something like that.)
    Extent <-250
    GrowWorld = matrix(0, nrow = Extent, ncol = Extent, byrow=TRUE)
    WeightWorld = matrix(0, nrow = Extent, ncol = Extent, byrow=TRUE)

    StartingSeed = as.data.frame(t(StartingSeed) )
    StartingSeed$row

    Recieved <- AddMoreSeedPoints(StartingSeed$row,StartingSeed$col, Patches, GrowWorld, WeightWorld, Extent)

    Patches <- Recieved$World
    GrowWorld <- Recieved$GrowWorld
    WeightWorld <- Recieved$WeightWorld
    #additional <- ifelse(Patches==max(Patches),1,0)
    #plot(raster(additional))
    World = Patches
    PatchesAndGrow = FillZeroPatch(Patches, GrowWorld, WeightWorld)
    Patches = PatchesAndGrow$Patches
    Possiblities = which(Patches==0 , arr.ind=TRUE)
  }
  return(Patches)
}



#' FillTheMatix adds habitats 2 to the number of MatrixHabs to the space between
#' habitat patches.
#'
#' This function adds habitats into the intervening space between a focal habitat patch.
#'
#' @param Patches array -
#' @param InterestLCV array -
#' @param MatrixHabs integer - default 10 habitats
#' @param PossibleNumberMatrixsPoints list of integers -
#' @param EqualSize boolean -
#' @param ChosenPatchDistribution character - how seed points are distributed to grow the patches. Options are "Random" or "Poisson", "Clustered" or "Thomas", or "Regular", "Strauss" or"Dispersed"
#' and defaults to random.
#' @param SeedDistribution character -a distribution from which to pick the size of patches from. Opetions are
#' "Gaussian" or"Normal", "Pareto", or the default "Uniform". Uniform gives patches of almost any size and the maximum variation.
#' @return list of arrays - the list contains two arrays, the "Patches" and "InterestLCV" (land cover).
#' @examples
#'
#' TBA...
#'
#' @export
FillTheMatix <- function(Patches, InterestLCV, MatrixHabs = 10, PossibleNumberMatrixsPoints = 1:200, EqualSize =FALSE, ChosenPatchDistribution="Random",SeedDistribution = "Uniform" ){

  if(length(PossibleNumberMatrixsPoints)==1){
    NumberOfEach = rep(PossibleNumberMatrixsPoints, MatrixHabs)
  }else{
    NumberOfEach = sample(PossibleNumberMatrixsPoints, MatrixHabs, replace = TRUE)
  }

  #print(NumberOfEach)
  UsedCells <- which(Patches>0, arr.ind=TRUE)
  MatrixCells = which(Patches==0, arr.ind=TRUE)
  MatrixSize = nrow(MatrixCells)

  if(EqualSize==TRUE){
    EachPatch = MatrixSize/sum(NumberOfEach)
    RandomEach = rep(EachPatch, sum(NumberOfEach))
  }else{
    RandomEach = sample(1:MatrixSize, sum(NumberOfEach), replace = TRUE)
  }


  Proportions = c()
  MaxSum = 0
  for(i in 1:length(NumberOfEach)){
    if(i ==1){
      MaxSum = NumberOfEach[i]
      Slicer = 1:MaxSum

    }else{
      Slicer = (1 + MaxSum):(MaxSum+NumberOfEach[i])
      MaxSum= MaxSum + NumberOfEach[i]
    }
    Proportions = c(Proportions,sum(RandomEach[Slicer]))
  }
  Sumproportions = sum(Proportions)
  Fract = Proportions/Sumproportions
  SizeEach = MatrixSize*Fract
  SizeEach = round(SizeEach)
  #print(SizeEach)
  MatrixHabitats = sample(2:11,MatrixHabs,replace = FALSE)

  Extent<- nrow(Patches)
  StartingCounter = length(Patches[Patches>0])
  EndingCounter = Extent^2
  NoInList = 1
  pb <- txtProgressBar(min = StartingCounter, max = EndingCounter, style = 3)
  for(hab in MatrixHabitats){
    setTxtProgressBar(pb, length(Patches[Patches > 0]))
    if (NoInList < length(MatrixHabitats)){
      #print(plot(raster(InterestLCV)))
      PatchAndLCV = CreateAdditionalWorld(Patches, InterestLCV, SizeEach[NoInList], NumberOfEach[NoInList], hab, PatchDistribution =ChosenPatchDistribution ,SeedDistribution)
      Patches = PatchAndLCV$Patches
      #print(plot(raster(Patches)))
      InterestLCV = PatchAndLCV$LCV
      #print(plot(raster(InterestLCV)))
    }else{
      Patches = FillAllZeroPatches(Patches)
      InterestLCV =ifelse(InterestLCV==0,hab, InterestLCV)
    }
    setTxtProgressBar(pb, length(Patches[Patches > 0]))
    NoInList = NoInList + 1
  }

  OutPut <- list("Patches" = Patches, "InterestLCV" = InterestLCV)
  return(OutPut)
}

#' FillTheMatix_specific adds habitats of specified type to the space between
#' habitat patches.
#'
#' This function adds habitats into the intervening space between a focal habitat patch.
#' It differs from \link[LcvGen]{FillTheMatix} by allowing non sequential habitats. For example,
#' habitat 2, 5, 7, 8, 9, 100.
#'
#' @param Patches array -
#' @param InterestLCV array -
#' @param MatrixHabs integer - default 10 habitats
#' @param PossibleNumberMatrixsPoints list of integers -
#' @param EqualSize boolean -
#' @param ChosenPatchDistribution character - how seed points are distributed to grow the patches. Options are "Random" or "Poisson", "Clustered" or "Thomas", or "Regular", "Strauss" or"Dispersed"
#' and defaults to random.
#' @param SeedDistribution character -a distribution from which to pick the size of patches from. Opetions are
#' "Gaussian" or"Normal", "Pareto", or the default "Uniform". Uniform gives patches of almost any size and the maximum variation.
#' @return list of arrays - the list contains two arrays, the "Patches" and "InterestLCV" (land cover).
#' @examples
#'
#' TBA...
#'
#' @export
FillTheMatix_specific <- function(Patches, InterestLCV, PossMatrixHabs = seq(2,11), PossibleNumberMatrixsPoints = 1:200, EqualSize =FALSE, ChosenPatchDistribution="Random",SeedDistribution = "Uniform" ){


  MatrixHabs <- length(PossMatrixHabs)



  if(length(PossibleNumberMatrixsPoints)==1){
    NumberOfEach = rep(PossibleNumberMatrixsPoints, MatrixHabs)
  }else{
    NumberOfEach = sample(PossibleNumberMatrixsPoints, MatrixHabs, replace = TRUE)
  }

  #print(NumberOfEach)
  UsedCells <- which(Patches>0, arr.ind=TRUE)
  MatrixCells = which(Patches==0, arr.ind=TRUE)
  MatrixSize = nrow(MatrixCells)

  if(EqualSize==TRUE){
    EachPatch = MatrixSize/sum(NumberOfEach)
    RandomEach = rep(EachPatch, sum(NumberOfEach))
  }else{
    RandomEach = sample(1:MatrixSize, sum(NumberOfEach), replace = TRUE)
  }


  Proportions = c()
  MaxSum = 0
  for(i in 1:length(NumberOfEach)){
    if(i ==1){
      MaxSum = NumberOfEach[i]
      Slicer = 1:MaxSum

    }else{
      Slicer = (1 + MaxSum):(MaxSum+NumberOfEach[i])
      MaxSum= MaxSum + NumberOfEach[i]
    }
    Proportions = c(Proportions,sum(RandomEach[Slicer]))
  }
  Sumproportions = sum(Proportions)
  Fract = Proportions/Sumproportions
  SizeEach = MatrixSize*Fract
  SizeEach = round(SizeEach)
  #print(SizeEach)
  MatrixHabitats = sample(PossMatrixHabs,MatrixHabs,replace = FALSE)

  Extent<- nrow(Patches)
  StartingCounter = length(Patches[Patches>0])
  EndingCounter = Extent^2
  NoInList = 1
  pb <- txtProgressBar(min = StartingCounter, max = EndingCounter, style = 3)
  for(hab in MatrixHabitats){
    setTxtProgressBar(pb, length(Patches[Patches > 0]))
    if (NoInList < length(MatrixHabitats)){
      #print(plot(raster(InterestLCV)))
      PatchAndLCV = CreateAdditionalWorld(Patches, InterestLCV, SizeEach[NoInList], NumberOfEach[NoInList], hab, PatchDistribution =ChosenPatchDistribution ,SeedDistribution)
      Patches = PatchAndLCV$Patches
      #print(plot(raster(Patches)))
      InterestLCV = PatchAndLCV$LCV
      #print(plot(raster(InterestLCV)))
    }else{
      Patches = FillAllZeroPatches(Patches)
      InterestLCV =ifelse(InterestLCV==0,hab, InterestLCV)
    }
    setTxtProgressBar(pb, length(Patches[Patches > 0]))
    NoInList = NoInList + 1
  }

  OutPut <- list("Patches" = Patches, "InterestLCV" = InterestLCV)
  return(OutPut)
}



ScaleUpArray <- function(InArray, scaleY, scaleX){
  OutArray = array(0, dim=c(nrow(InArray)*scaleY, ncol(InArray)*scaleX))
  for(y in 1:nrow(InArray)){
    for(x in 1:ncol(InArray)){
      for(iy in (y*scaleY-(scaleY-1)):(y*scaleY)){
        for(ix in (x*scaleX-(scaleX-1)):(x*scaleX)){
          OutArray[iy, ix]=InArray[y,x]
        }
      }
    }
  }
  return(OutArray)
}

AddBuffer = function(landcovermatrix, patchesmatrix, buffer=5){


  newrows = nrow(landcovermatrix) + (2 * buffer)
  newcols = ncol(landcovermatrix) + (2 * buffer)

  myMatrix = matrix(data=0, nrow = newrows, ncol = newcols, byrow = TRUE)
  OutLandcovermatrix = myMatrix

  OutPatchesmatrix = myMatrix



  for (i in 1:nrow(landcovermatrix)){
    for( j in 1:ncol(landcovermatrix)){
      OutLandcovermatrix[j+buffer, i+buffer]= landcovermatrix[j,i]
      OutPatchesmatrix[j+buffer, i+buffer]= patchesmatrix[j,i]
    }
  }

  OutBoth = list("OutPatchesmatrix" = OutPatchesmatrix, "OutLandcovermatrix"=OutLandcovermatrix)

}




#' GenerateLandcover creates a ascii land cover dataset with a corresponding patch number data set.
#'
#' This function uses other functions to create a land cover with the focal habitat haveing a number
#' of patches.
#'
#' @param NumberSeeds integer -
#' @param LCVPath character - string of the output path of the land cover ascii
#' @param PatchPath character - string of the output path of the patch ascii
#' @param runno integer - defaults to NULL
#' @param EqualSize boolean -
#' @param ChosenPatchDistribution character - how seed points are distributed to grow the patches. Options are "Random" or "Poisson", "Clustered" or "Thomas", or "Regular", "Strauss" or"Dispersed"
#' and defaults to random.
#' @param SeedDistribution character -a distribution from which to pick the size of patches from. Opetions are
#' "Gaussian" or"Normal", "Pareto", or the default "Uniform". Uniform gives patches of almost any size and the maximum variation.
#' @return list of arrays - the list contains two arrays, the "Patches" and "InterestLCV" (land cover).
#' @examples
#'
#' TBA...
#'
#' @export
GenerateLandcover <-function(NumberSeeds, LCVPath, PatchPath, runno = NULL, MatrixHabs =10){
  #require(raster)
  #require(rasterVis)
  "
  This sets up the habitat of interest patches.

  We can very how many patches.

  They can be roughly the same size, or their sizes can be random (probability of growth from a uniform
  distribution).

  The starting points can be:
  #Clustered
  #Random
  #Distributed

  Adittionally two different percentage covers can be generated from the same seed locations.

  "

  TheExtent = 250
  #axisLabels = list(at=c(0, 250, 500, 750,1000))

  InterestPatches = CreateWorld(SeedPoints = NumberSeeds, PercentageCover = 0.1, Extent = TheExtent, SecondPercentageCover = 0 , PatchDistribution = "Random", SeedDistribution = "Uniform")

  Patches = InterestPatches$LargerWorld
  LCV = ifelse(Patches>0,1,0)

  PatchesDispersedMatrix <- FillTheMatix(Patches, LCV, MatrixHabs = MatrixHabs, PossibleNumberMatrixsPoints = 1:200, EqualSize = FALSE, ChosenPatchDistribution = "Random", SeedDistribution = "Uniform")

  Patches = PatchesDispersedMatrix$Patches
  LCV = PatchesDispersedMatrix$InterestLCV

  Patches = ScaleUpArray(Patches,4,4)
  LCV = ScaleUpArray(LCV,4,4)

  WithBuffer = AddBuffer(LCV, Patches, 10)
  Patches = WithBuffer$OutPatchesmatrix
  LCV = WithBuffer$OutLandcovermatrix

  #RasPatches = raster(Patches, xmn=0, xmx =1000, ymn=0, ymx=1000)
  #levelplot(RasPatches , col.regions = rev(terrain.colors(255)), margin=FALSE, scales=list(x=axisLabels,y=axisLabels))

  #RasLCV = raster(LCV, xmn=0, xmx =1000, ymn=0, ymx=1000)
  #levelplot(RasLCV, col.regions = rev(terrain.colors(255)), margin=FALSE, scales=list(x=axisLabels,y=axisLabels))

  ArrayToAsc(Patches, filename=PatchPath, overwrite = TRUE, runno)
  #writeRaster(RasPatches, filename=PatchPath, datatype="ascii", overwrite=TRUE)
  #writeRaster( RasLCV, filename=LCVPath, datatype="ascii", overwrite=TRUE)
  ArrayToAsc(LCV, filename=LCVPath, overwrite = TRUE, runno)
  OutList = list("Patches"=Patches, "LCV"=LCV )
  return(OutList)
}

#'To make this a package it needs to also have the capability to specify possible habitats and if they are all used, or how many are.'

#Out = GenerateLandcover(10, "LCV_test2.asc", "Patch_test2.asc", runno = 2)

# parameters <- read.csv('parameters.csv', stringsAsFactors=FALSE)
#
# for(i in 1:10400){
#   lcver = paste("E:/SimedRasters/LCV", i, ".asc")
#   patches = paste("E:/SimedRasters/PATCHES", i, ".asc")
#
#   NoPatches = parameters[which(parameters$RunNo==i),]$NoPatches
#   GenerateLandcover(NoPatches, lcver, patches)
# }

