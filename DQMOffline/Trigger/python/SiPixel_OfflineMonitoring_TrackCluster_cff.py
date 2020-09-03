import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester

from DQMOffline.Trigger.SiPixel_OfflineMonitoring_HistogramManager_cfi import *

# order is important and it should follow ordering in hltSiPixelPhase1ClustersConf VPSet
hltSiPixelPhase1TrackClustersOnTrackCharge = hltDefaultHistoTrack.clone(
  name = "charge",
  title = "Corrected Cluster Charge (OnTrack)",
  range_min = 0, range_max = 200e3, range_nbins = 200,
  xlabel = "Charge (electrons)",
  specs = VPSet(
    hltStandardSpecifications1D,
    StandardSpecification2DProfile
  )
)

hltSiPixelPhase1TrackClustersOffTrackCharge = hltDefaultHistoTrack.clone(
  name = "charge_OffTrack",
  title = "Corrected Cluster Charge (OffTrack)",
  range_min = 0, range_max = 200e3, range_nbins = 200,
  xlabel = "Charge (electrons)",
  specs = VPSet(
    hltStandardSpecifications1D,
    StandardSpecification2DProfile
  )
)

hltSiPixelPhase1TrackClustersOnTrackBigPixelCharge = hltDefaultHistoTrack.clone(
  name = "bigpixelcharge",
  title = "Corrected Big Pixel Charge (OnTrack)",
  range_min = 0, range_max = 80e3, range_nbins = 100,
  xlabel = "Charge (electrons)",
  enabled=False,

  specs = VPSet(
    Specification().groupBy("PXBarrel").save(),
    Specification().groupBy("PXForward").save(),
    Specification().groupBy("PXBarrel/PXLayer").save(),
    Specification().groupBy("PXForward/PXDisk").save()
  )
)

hltSiPixelPhase1TrackClustersOffTrackBigPixelCharge = hltDefaultHistoTrack.clone(
  name = "bigpixelcharge_OffTrack",
  title = "Corrected Big Pixel Charge (OffTrack)",
  range_min = 0, range_max = 80e3, range_nbins = 100,
  xlabel = "Charge (electrons)",
  enabled=False,

  specs = VPSet(
    Specification().groupBy("PXBarrel").save(),
    Specification().groupBy("PXForward").save(),
    Specification().groupBy("PXBarrel/PXLayer").save(),
    Specification().groupBy("PXForward/PXDisk").save()
  )
)

hltSiPixelPhase1TrackClustersOnTrackNotBigPixelCharge = hltDefaultHistoTrack.clone(
  name = "notbigpixelcharge",
  title = "Corrected Not Big Pixel Charge (OnTrack)",
  range_min = 0, range_max = 80e3, range_nbins = 100,
  xlabel = "Charge (electrons)",
  enabled=False,

  specs = VPSet(
    Specification().groupBy("PXBarrel").save(),
    Specification().groupBy("PXForward").save(),
    Specification().groupBy("PXBarrel/PXLayer").save(),
    Specification().groupBy("PXForward/PXDisk").save()
  )
)

hltSiPixelPhase1TrackClustersOffTrackNotBigPixelCharge = hltDefaultHistoTrack.clone(
  name = "notbigpixelcharge_OffTrack",
  title = "Corrected Not Big Pixel Charge (OffTrack)",
  range_min = 0, range_max = 80e3, range_nbins = 100,
  xlabel = "Charge (electrons)",
  enabled=False,

  specs = VPSet(
    Specification().groupBy("PXBarrel").save(),
    Specification().groupBy("PXForward").save(),
    Specification().groupBy("PXBarrel/PXLayer").save(),
    Specification().groupBy("PXForward/PXDisk").save()
  )
)

hltSiPixelPhase1TrackClustersOnTrackSize = hltDefaultHistoTrack.clone(
  name = "size",
  title = "Total Cluster Size (OnTrack)",
  range_min = 0, range_max = 30, range_nbins = 30,
  xlabel = "size[pixels]",

  specs = VPSet(
    hltStandardSpecifications1D
  )
)

hltSiPixelPhase1TrackClustersOffTrackSize = hltDefaultHistoTrack.clone(
  name = "size_OffTrack",
  title = "Total Cluster Size (OffTrack)",
  range_min = 0, range_max = 30, range_nbins = 30,
  xlabel = "size[pixels]",

  specs = VPSet(
    hltStandardSpecifications1D
  )
)

hltSiPixelPhase1TrackClustersOffTrackSizex = hltDefaultHistoTrack.clone(
  name = "sizex_OffTrack",
  title = "Cluster Size along radial (OffTrack)",
  range_min = 0, range_max = 30, range_nbins = 30,
  xlabel = "size[pixels]",

  specs = VPSet(
    hltStandardSpecifications1D
  )
)

hltSiPixelPhase1TrackClustersOffTrackSizey = hltDefaultHistoTrack.clone(
  name = "sizey_OffTrack",
  title = "Cluster Size along Beamline (OffTrack)",
  range_min = 0, range_max = 30, range_nbins = 30,
  xlabel = "size[pixels]",

  specs = VPSet(
    hltStandardSpecifications1D
  )
)

hltSiPixelPhase1TrackClustersOnTrackShape = hltDefaultHistoTrack.clone(
  enabled = False,
  name = "shapeFilter",
  title = "Shape filter (OnTrack)",
  range_min = 0, range_max = 2, range_nbins = 2,
  xlabel = "shapeFilter",

  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer").saveAll(),
    Specification().groupBy("PXForward/PXDisk").saveAll(),
    StandardSpecification2DProfile,

    Specification().groupBy("PXBarrel/PXLayer/Lumisection")
                   .reduce("MEAN")
                   .groupBy("PXBarrel/PXLayer", "EXTEND_X")
                   .save(),

    Specification().groupBy("PXForward/PXDisk/Lumisection")
                   .reduce("MEAN")
                   .groupBy("PXForward/PXDisk", "EXTEND_X")
                   .save(),
  )
)

hltSiPixelPhase1TrackClustersOnTrackNClusters = hltDefaultHistoTrack.clone(
  name = "clusters_ontrack",
  title = "Clusters_onTrack",
  range_min = 0, range_max = 40000, range_nbins = 4000,
  xlabel = "clusters",
  dimensions = 0,
  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer" + "/DetId/Event")
                   .reduce("COUNT")
                   .groupBy("PXBarrel/PXLayer")
                   .saveAll(),
    Specification().groupBy("PXForward/PXDisk" + "/DetId/Event")
                   .reduce("COUNT")
                   .groupBy("PXForward/PXDisk")
                   .saveAll(),
    StandardSpecificationInclusive_Num,
    StandardSpecificationTrend_Num
  )
)

hltSiPixelPhase1TrackClustersOffTrackNClusters = hltDefaultHistoTrack.clone(
  name = "clusters_OffTrack",
  title = "Clusters_offTrack",
  range_min = 0, range_max = 50000, range_nbins = 5000,
  xlabel = "clusters",
  dimensions = 0,
  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer" + "/DetId/Event")
                   .reduce("COUNT")
                   .groupBy("PXBarrel/PXLayer")
                   .saveAll(),
    Specification().groupBy("PXForward/PXDisk" + "/DetId/Event")
                   .reduce("COUNT")
                   .groupBy("PXForward/PXDisk")
                   .saveAll(),
    StandardSpecificationInclusive_Num,
    StandardSpecificationTrend_Num
  )
)

hltSiPixelPhase1TrackClustersOnTrackPositionB = hltDefaultHistoTrack.clone(
  name = "clusterposition_zphi_ontrack",
  title = "Cluster_onTrack Positions",
  range_min   =  -60, range_max   =  60, range_nbins   = 600,
  range_y_min = -3.2, range_y_max = 3.2, range_y_nbins = 200,
  xlabel = "Global Z", ylabel = "Global \phi",
  dimensions = 2,
  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer").save(),
#    Specification().groupBy("").save(),
  )
)

hltSiPixelPhase1TrackClustersOffTrackPositionB = hltDefaultHistoTrack.clone(
  name = "clusterposition_zphi_OffTrack",
  title = "Cluster_offTrack Positions",
  range_min   =  -60, range_max   =  60, range_nbins   = 600,
  range_y_min = -3.2, range_y_max = 3.2, range_y_nbins = 200,
  xlabel = "Global Z", ylabel = "Global \phi",
  dimensions = 2,
  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer").save(),
#    Specification().groupBy("").save(),
  )
)

hltSiPixelPhase1TrackClustersOnTrackPositionF = hltDefaultHistoTrack.clone(
  name = "clusterposition_xy_ontrack",
  title = "Cluster_onTrack Positions",
  xlabel = "Global X", ylabel = "Global Y",
  range_min   = -20, range_max   = 20, range_nbins   = 200,
  range_y_min = -20, range_y_max = 20, range_y_nbins = 200,
  dimensions = 2,
  specs = VPSet(
    Specification().groupBy("PXForward/PXDisk").save(),
  )
)

hltSiPixelPhase1TrackClustersOffTrackPositionF = hltDefaultHistoTrack.clone(
  name = "clusterposition_xy_offtrack",
  title = "Cluster_offTrack Positions",
  xlabel = "Global X", ylabel = "Global Y",
  range_min   = -20, range_max   = 20, range_nbins   = 200,
  range_y_min = -20, range_y_max = 20, range_y_nbins = 200,
  dimensions = 2,
  specs = VPSet(
    Specification().groupBy("PXForward/PXDisk").save(),
  )
)

hltSiPixelPhase1TrackClustersOffTrackPositionXz = hltDefaultHistoTrack.clone(
  name = "clusterposition_xz_OffTrack",
  title = "Cluster_offTrack Positions",
  xlabel = "Global X", ylabel = "Global Z",
  range_min   = -20, range_max   = 20, range_nbins   = 200,
  range_y_min = -20, range_y_max = 20, range_y_nbins = 200,
  dimensions = 2,
  specs = VPSet(
    Specification().groupBy("PXForward/PXDisk").save(),
  )
)

hltSiPixelPhase1TrackClustersOffTrackPositionYz = hltDefaultHistoTrack.clone(
  name = "clusterposition_yz_OffTrack",
  title = "Cluster_offTrack Positions",
  xlabel = "Global y", ylabel = "Global Z",
  range_min   = -20, range_max   = 20, range_nbins   = 200,
  range_y_min = -20, range_y_max = 20, range_y_nbins = 200,
  dimensions = 2,
  specs = VPSet(
    Specification().groupBy("PXForward/PXDisk").save(),
  )
)

hltSiPixelPhase1DigisHitmapOnTrack = hltDefaultHistoTrack.clone(
  enabled = False,
  name = "digi_occupancy_ontrack",
  title = "Digi Occupancy (OnTrack)",
  ylabel = "#digis",
  dimensions = 0,
  specs = VPSet(
    Specification(PerModule).groupBy("PXBarrel/Shell/PXLayer/SignedLadder/PXModuleName/row/col")
                            .groupBy("PXBarrel/Shell/PXLayer/SignedLadder/PXModuleName/row", "EXTEND_X")
                            .groupBy("PXBarrel/Shell/PXLayer/SignedLadder/PXModuleName", "EXTEND_Y")
                            .save(),
    Specification(PerModule).groupBy("PXForward/HalfCylinder/PXRing/PXDisk/SignedBlade/PXModuleName/row/col")
                            .groupBy("PXForward/HalfCylinder/PXRing/PXDisk/SignedBlade/PXModuleName/row", "EXTEND_X")
                            .groupBy("PXForward/HalfCylinder/PXRing/PXDisk/SignedBlade/PXModuleName", "EXTEND_Y")
                            .save(),
  )
)

hltSiPixelPhase1DigisHitmapOffTrack = hltDefaultHistoTrack.clone(
  enabled = False,
  name = "digi_occupancy_OffTrack",
  title = "Digi Occupancy (OffTrack)",
  ylabel = "#digis",
  dimensions = 0,
  specs = VPSet(
    Specification(PerModule).groupBy("PXBarrel/Shell/PXLayer/SignedLadder/PXModuleName/row/col")
                            .groupBy("PXBarrel/Shell/PXLayer/SignedLadder/PXModuleName/row", "EXTEND_X")
                            .groupBy("PXBarrel/Shell/PXLayer/SignedLadder/PXModuleName", "EXTEND_Y")
                            .save(),
    Specification(PerModule).groupBy("PXForward/HalfCylinder/PXRing/PXDisk/SignedBlade/PXModuleName/row/col")
                            .groupBy("PXForward/HalfCylinder/PXRing/PXDisk/SignedBlade/PXModuleName/row", "EXTEND_X")
                            .groupBy("PXForward/HalfCylinder/PXRing/PXDisk/SignedBlade/PXModuleName", "EXTEND_Y")
                            .save(),
  )
)

hltSiPixelPhase1DigisNdigisOnTrack = hltDefaultHistoTrack.clone(
  enabled = False,
  name = "digis on-track", # 'Count of' added automatically
  title = "Digis on-track",
  xlabel = "digis (on-track)",
  range_min = 0,
  range_max = 300,
  range_nbins = 50,
  dimensions = 0, # this is a count

  specs = VPSet(
    StandardSpecificationTrend_Num,
    StandardSpecification2DProfile_Num,
    
    Specification().groupBy("PXBarrel/PXLayer/Event") #this will produce inclusive counts per Layer/Disk
                             .reduce("COUNT")    
                             .groupBy("PXBarrel/PXLayer")
                             .save(nbins=100, xmin=0, xmax=40000),
    Specification().groupBy("PXForward/PXDisk/Event")
                             .reduce("COUNT")    
                             .groupBy("PXForward/PXDisk/")
                             .save(nbins=100, xmin=0, xmax=20000),
  )
)

hltSiPixelPhase1DigisNdigisOffTrack = hltDefaultHistoTrack.clone(
  enabled = False,
  name = "digis_OffTrack",# 'Count of' added automatically
  title = "Digis off-track",
  xlabel = "digis (off-track)",
  range_min = 0,
  range_max = 300,
  range_nbins = 50,
  dimensions = 0, # this is a count

  specs = VPSet(
    StandardSpecificationTrend_Num,
    StandardSpecification2DProfile_Num,
    
    Specification().groupBy("PXBarrel/PXLayer/Event") #this will produce inclusive counts per Layer/Disk
                             .reduce("COUNT")    
                             .groupBy("PXBarrel/PXLayer")
                             .save(nbins=100, xmin=0, xmax=40000),
    Specification().groupBy("PXForward/PXDisk/Event")
                             .reduce("COUNT")    
                             .groupBy("PXForward/PXDisk/")
                             .save(nbins=100, xmin=0, xmax=20000),
  )
)


hltSiPixelPhase1TrackClustersDigisOverClusterTotcharge = hltDefaultHistoTrack.clone(
  name = "digis_cluster_chargeratio",
  title = "TotalDigis/Cluster Charge ratio Digis",
  xlabel = "DigisTotalCharge/ClusterCharge",
  range_min = 0,
  range_max = 100,
  range_nbins = 100,
  dimensions = 1,

  specs = VPSet(
    StandardSpecifications1D,
    #Specification().groupBy("PXBarrel/PXLayer").save(),
    #Specification().groupBy("PXForward/PXDisk").save(),
  )
)

hltSiPixelPhase1TrackClustersDigisOverClusterTotcharge2D = hltDefaultHistoTrack.clone(
  name = "total_digis_charge_vs_cluster_charge",
  title = "TotalDigis Charge ratio vs Cluster Charge",
  ylabel = "Digis Total Charge",
  xlabel = "Cluster Charge",
  range_y_min = 0,
  range_y_max = 100,
  range_y_nbins = 100,
  range_min = 0, range_max = 80e3, range_nbins = 100,
  dimensions = 2,

  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer").save(),
    Specification().groupBy("PXForward/PXDisk").save()
  )
)

hltSiPixelPhase1TrackClustersNTracks = hltDefaultHistoTrack.clone(
  name = "ntracks",
  title = "Number of Tracks",
  xlabel = "All - Pixel - BPIX - FPIX",
  range_min = 1, range_max = 5, range_nbins = 4,
  dimensions = 1,
  specs = VPSet(
    Specification().groupBy("").save()
  )
)

hltSiPixelPhase1TrackClustersNTracksInVolume = hltDefaultHistoTrack.clone(
  name = "ntracksinpixvolume",
  title = "Number of Tracks in Pixel fiducial Volume",
  xlabel = "without hits - with hits",
  range_min = 0, range_max = 2, range_nbins = 2,
  dimensions = 1,
  specs = VPSet(
    Specification().groupBy("").save()
  )
)

hltSiPixelPhase1ClustersSizeVsEtaOnTrackOuter = hltDefaultHistoTrack.clone(
  name = "sizeyvseta_on_track_outer",
  title = "Cluster Size along Beamline vs. Track #eta (OnTrack) outer ladders",
  xlabel = "Track #eta",
  ylabel = "length [pixels]",
  range_min = -3.2, range_max  = 3.2, range_nbins   = 64,
  range_y_min =  0, range_y_max = 30, range_y_nbins = 30,
  dimensions = 2,
  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer").save()
  )
)

hltSiPixelPhase1ClustersSizeVsEtaOnTrackInner = hltSiPixelPhase1ClustersSizeVsEtaOnTrackOuter.clone(
  name = "sizeyvseta_on_track_inner",
  title = "Cluster Size along Beamline vs. Track #eta (OnTrack) inner ladders",
)


hltSiPixelPhase1TrackClustersOnTrackSizeYOuter = hltSiPixelPhase1ClustersSizeVsEtaOnTrackOuter.clone(
  name = "sizey_on_track_outer",
  title = "Cluster Size along Beamline vs. prediction (OnTrack) outer ladders",
  xlabel = "prediction",
  ylabel = "length [pixels]",
  range_min = 0, range_max  = 30, range_nbins   = 60
)

hltSiPixelPhase1TrackClustersOnTrackSizeYInner = hltSiPixelPhase1TrackClustersOnTrackSizeYOuter.clone(
  name = "sizey_on_track_inner",
  title = "Cluster Size along Beamline vs. prediction (OnTrack) inner ladders",
)

hltSiPixelPhase1TrackClustersOnTrackSizeYF = hltSiPixelPhase1TrackClustersOnTrackSizeYOuter.clone(
  name = "sizey_on_track_forward",
  title = "Cluster Size ortogonal to Beamline vs. prediction (OnTrack) forward",
  range_y_min =  0, range_y_max = 10, range_y_nbins = 10,
  range_min = 0, range_max  = 10, range_nbins   = 20,
  specs = VPSet(
    Specification().groupBy("PXForward/PXDisk").save(),
  )
)


hltSiPixelPhase1TrackClustersOnTrackSizeXOuter = hltSiPixelPhase1TrackClustersOnTrackSizeYOuter.clone(
  name = "sizex_on_track_outer",
  title = "Cluster Size along radial vs. prediction (OnTrack) outer ladders",
  range_min = 0, range_max  = 8, range_nbins   = 16,
  range_y_min =  0, range_y_max = 8, range_y_nbins = 8

)

hltSiPixelPhase1TrackClustersOnTrackSizeXInner = hltSiPixelPhase1TrackClustersOnTrackSizeXOuter.clone(
  name = "sizex_on_track_inner",
  title = "Cluster Size along radial vs. prediction (OnTrack) inner ladders",
)

hltSiPixelPhase1TrackClustersOnTrackSizeXF = hltSiPixelPhase1TrackClustersOnTrackSizeYF.clone(
  name = "sizex_on_track_forward",
  title = "Cluster Size radial vs. prediction (OnTrack) forward",
)



hltSiPixelPhase1TrackClustersOnTrackSizeXYOuter = hltSiPixelPhase1TrackClustersOnTrackSizeYOuter.clone(
  name = "sizexy_on_track_outer",
  title = "Cluster Size x vs y (OnTrack) outer ladders",
  xlabel = "y size",
  ylabel = "x size",
  range_min = 0, range_max  = 20, range_nbins   = 20,
  range_y_min = 0, range_y_max = 10, range_y_nbins = 10
)

hltSiPixelPhase1TrackClustersOnTrackSizeXYInner = hltSiPixelPhase1TrackClustersOnTrackSizeXYOuter.clone(
 name = "sizexy_on_track_inner",
 title = "Cluster Size x vs y (OnTrack) inner ladders"
)

hltSiPixelPhase1TrackClustersOnTrackSizeXYF = hltSiPixelPhase1TrackClustersOnTrackSizeYF.clone(
  name = "sizexy_on_track_forward",
  title = "Cluster Size x vs y (OnTrack) forward",
  xlabel = "y size",
  ylabel = "x size",
  range_min = 0, range_max  = 10, range_nbins   = 10,
  range_y_min = 0, range_y_max = 10, range_y_nbins = 10

)

hltSiPixelPhase1ClustersChargeVsSizeOnTrack = hltDefaultHistoTrack.clone(
  name = "chargevssize_on_track",
  title = "Cluster Charge vs. Cluster Size (OnTrack)",
  xlabel = "size[pixels]",
  ylabel = "Cluster charge",
  range_min =  0, range_max = 30, range_nbins = 15,
  range_y_min = 0, range_y_max = 80e3, range_y_nbins = 100,
  dimensions = 2,
  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer").save(),
    Specification().groupBy("PXForward/PXDisk").save()
  )
)

hltSiPixelPhase1ClustersSizeVsEtaOffTrack = hltDefaultHistoTrack.clone(
  name = "sizeyvseta__OffTrack",
  title = "Cluster Size along Beamline vs. Cluster #eta (OffTrack)",
  xlabel = "Cluster #eta",
  ylabel = "length [pixels]",
  range_min = -3.2, range_max  = 3.2, range_nbins   = 64,
  range_y_min =  0, range_y_max = 30, range_y_nbins = 15,
  dimensions = 2,
  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer").save(),
    Specification().groupBy("PXForward/PXDisk").save()
  )
)

hltSiPixelPhase1ClustersChargeVsEtaOffTrack = DefaultHistoTrack.clone(
  name = "Chargevseta_OffTrack",
  title = "Cluster Charge vs. Cluster #eta (OffTrack)",
  xlabel = "Cluster #eta",
  ylabel = "Cluster Charge",
  range_min = -3.2, range_max  = 3.2, range_nbins   = 64,
  range_y_min = 0, range_y_max = 80e3, range_y_nbins = 100,
  dimensions = 2,
  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer").save(),
    Specification().groupBy("PXForward/PXDisk").save()
  )
)

hltSiPixelPhase1ClustersChargeVsSizeOffTrack = hltDefaultHistoTrack.clone(
  name = "Chargevssizey_OffTrack",
  title = "Cluster Charge vs. Cluster Size along Beamline",
  xlabel = "length [pixels]",
  ylabel = "Cluster Charge",
  range_min =  0, range_max = 30, range_nbins = 15,
  range_y_min = 0, range_y_max = 80e3, range_y_nbins = 100,
  dimensions = 2,
  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer").save(),
    Specification().groupBy("PXForward/PXDisk").save()
  )
)

hltSiPixelPhase1TrackClustersOnTrackChargeOuter = hltDefaultHistoTrack.clone(
  name = "chargeOuter",
  title = "Corrected Cluster Charge (OnTrack) outer ladders",
  range_min = 0, range_max = 150e3, range_nbins = 150,
  xlabel = "Charge (electrons)",

  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer").save()
  )
)

hltSiPixelPhase1TrackClustersOnTrackChargeInner = hltSiPixelPhase1TrackClustersOnTrackChargeOuter.clone(
  name = "chargeInner",
  title = "Corrected Cluster Charge (OnTrack) inner ladders"
)

hltSiPixelPhase1TrackClustersOnTrackShapeOuter = hltDefaultHistoTrack.clone(
  enabled = False,
  name = "shapeFilterOuter",
  title = "Shape filter (OnTrack) Outer Ladders",
  range_min = 0, range_max = 2, range_nbins = 2,
  xlabel = "shapeFilter",
  specs = VPSet(
    Specification().groupBy("PXBarrel/PXLayer").save()
  )
)

hltSiPixelPhase1TrackClustersOnTrackShapeInner = hltSiPixelPhase1TrackClustersOnTrackShapeOuter.clone(
  enabled = False,
  name = "shapeFilterInner",
  title = "Shape filter (OnTrack) Inner Ladders",
)

hltSiPixelPhase1TrackClustersConf = cms.VPSet(
### IT HAS TO BE IN SYNCH W/
### THE LIST DEFINED IN THE ENUM
### https://cmssdt.cern.ch/lxr/source/DQM/SiPixelPhase1TrackClusters/src/SiPixelPhase1TrackClusters.cc#0063
   hltSiPixelPhase1TrackClustersOnTrackCharge,
   hltSiPixelPhase1TrackClustersOffTrackCharge,
   hltSiPixelPhase1TrackClustersOnTrackBigPixelCharge,
   hltSiPixelPhase1TrackClustersOffTrackBigPixelCharge,
   hltSiPixelPhase1TrackClustersOnTrackNotBigPixelCharge,
   hltSiPixelPhase1TrackClustersOffTrackNotBigPixelCharge,

   hltSiPixelPhase1TrackClustersOnTrackSize,
   hltSiPixelPhase1TrackClustersOffTrackSize,
   hltSiPixelPhase1TrackClustersOffTrackSizex,
   hltSiPixelPhase1TrackClustersOffTrackSizey,
   hltSiPixelPhase1TrackClustersOnTrackShape,

   hltSiPixelPhase1TrackClustersOnTrackNClusters,
   hltSiPixelPhase1TrackClustersOffTrackNClusters,

   hltSiPixelPhase1TrackClustersOnTrackPositionB,
   hltSiPixelPhase1TrackClustersOffTrackPositionB,
   hltSiPixelPhase1TrackClustersOnTrackPositionF,
   hltSiPixelPhase1TrackClustersOffTrackPositionF,
   hltSiPixelPhase1TrackClustersOffTrackPositionXz,
   hltSiPixelPhase1TrackClustersOffTrackPositionYz,

   hltSiPixelPhase1DigisHitmapOnTrack,
   hltSiPixelPhase1DigisHitmapOffTrack,
   hltSiPixelPhase1DigisNdigisOnTrack,
   hltSiPixelPhase1DigisNdigisOffTrack,

   hltSiPixelPhase1TrackClustersDigisOverClusterTotcharge,
   hltSiPixelPhase1TrackClustersDigisOverClusterTotcharge2D,

   hltSiPixelPhase1TrackClustersNTracks,
   hltSiPixelPhase1TrackClustersNTracksInVolume,

   hltSiPixelPhase1ClustersSizeVsEtaOnTrackOuter,
   hltSiPixelPhase1ClustersSizeVsEtaOnTrackInner,
   hltSiPixelPhase1TrackClustersOnTrackChargeOuter,
   hltSiPixelPhase1TrackClustersOnTrackChargeInner,

   hltSiPixelPhase1TrackClustersOnTrackShapeOuter,
   hltSiPixelPhase1TrackClustersOnTrackShapeInner,

   hltSiPixelPhase1TrackClustersOnTrackSizeXOuter,
   hltSiPixelPhase1TrackClustersOnTrackSizeXInner,
   hltSiPixelPhase1TrackClustersOnTrackSizeXF,
   hltSiPixelPhase1TrackClustersOnTrackSizeYOuter,
   hltSiPixelPhase1TrackClustersOnTrackSizeYInner,
   hltSiPixelPhase1TrackClustersOnTrackSizeYF,

   hltSiPixelPhase1TrackClustersOnTrackSizeXYOuter,
   hltSiPixelPhase1TrackClustersOnTrackSizeXYInner,
   hltSiPixelPhase1TrackClustersOnTrackSizeXYF,
   hltSiPixelPhase1ClustersChargeVsSizeOnTrack,
   hltSiPixelPhase1ClustersSizeVsEtaOffTrack,
   hltSiPixelPhase1ClustersChargeVsEtaOffTrack,
   hltSiPixelPhase1ClustersChargeVsSizeOffTrack,
)

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
hltSiPixelPhase1TrackClustersAnalyzer = DQMEDAnalyzer('SiPixelPhase1TrackClusters',
        VertexCut  = cms.untracked.bool(False),
        clusters   = cms.InputTag("hltSiPixelClusters"),
        tracks     = cms.InputTag("hltMergedTracks"), #hltIter2Merged"
        clusterShapeCache = cms.InputTag("hltSiPixelClusterShapeCache"),
        vertices = cms.InputTag(""),
        histograms = hltSiPixelPhase1TrackClustersConf,
        geometry   = hltSiPixelPhase1Geometry
)

hltSiPixelPhase1TrackClustersHarvester = DQMEDHarvester("SiPixelPhase1Harvester",
        histograms = hltSiPixelPhase1TrackClustersConf,
        geometry   = hltSiPixelPhase1Geometry
)
