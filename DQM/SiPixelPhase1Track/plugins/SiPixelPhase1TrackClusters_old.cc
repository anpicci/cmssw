// -*- C++ -*-
//
// Package:     SiPixelPhase1TrackClusters
// Class  :     SiPixelPhase1TrackClusters
//

// Original Author: Marcel Schneider

#include "DQM/SiPixelPhase1Common/interface/SiPixelPhase1Base.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelClusterShapeCache.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"

#include "RecoTracker/Record/interface/CkfComponentsRecord.h"
#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterShapeHitFilter.h"

#include <iostream>

namespace {

  class SiPixelPhase1TrackClusters final : public SiPixelPhase1Base {
    enum {
      ON_TRACK_CHARGE,
      OFF_TRACK_CHARGE,
      ON_TRACK_BIGPIXELCHARGE,
      OFF_TRACK_BIGPIXELCHARGE,
      ON_TRACK_NOTBIGPIXELCHARGE,
      OFF_TRACK_NOTBIGPIXELCHARGE,

      ON_TRACK_SIZE,
      OFF_TRACK_SIZE,
      OFF_TRACK_SIZEX,
      OFF_TRACK_SIZEY,
      ON_TRACK_SHAPE,

      ON_TRACK_NCLUSTERS,
      OFF_TRACK_NCLUSTERS,

      ON_TRACK_POSITIONB,
      OFF_TRACK_POSITION_B,
      ON_TRACK_POSITIONF,
      OFF_TRACK_POSITION_F,
      OFF_TRACK_POSITION_XZ,
      OFF_TRACK_POSITION_YZ,

      DIGIS_HITMAP_ON_TRACK,
      DIGIS_HITMAP_OFF_TRACK,
      ON_TRACK_NDIGIS,
      OFF_TRACK_NDIGIS,

      //DIGIS_OVER_CLUSTER_TOTCHARGE,
      //DIGIS_OVER_CLUSTER_TOTCHARGE_2D,

      //OFF_TRACK_READOUT_CHARGE,
      //OFF_TRACK_READOUT_NCLUSTERS,
      
      NTRACKS,
      NTRACKS_INVOLUME,

      SIZE_VS_ETA_ON_TRACK_OUTER,
      SIZE_VS_ETA_ON_TRACK_INNER,
      ON_TRACK_CHARGE_OUTER,
      ON_TRACK_CHARGE_INNER,

      ON_TRACK_SHAPE_OUTER,
      ON_TRACK_SHAPE_INNER,

      ON_TRACK_SIZE_X_OUTER,
      ON_TRACK_SIZE_X_INNER,
      ON_TRACK_SIZE_X_F,
      ON_TRACK_SIZE_Y_OUTER,
      ON_TRACK_SIZE_Y_INNER,
      ON_TRACK_SIZE_Y_F,

      ON_TRACK_SIZE_XY_OUTER,
      ON_TRACK_SIZE_XY_INNER,
      ON_TRACK_SIZE_XY_F,
      CHARGE_VS_SIZE_ON_TRACK,

      SIZE_VS_ETA_OFF_TRACK,
      CHARGE_VS_ETA_OFF_TRACK,
      CHARGE_VS_SIZE_OFF_TRACK,
      
      ENUM_SIZE
    };

  public:
    explicit SiPixelPhase1TrackClusters(const edm::ParameterSet& conf);
    //explicit SiPixelPhase1Clusters(const edm::ParameterSet& conf);
    void analyze(const edm::Event&, const edm::EventSetup&) override;

  private:
    const bool applyVertexCut_;
    edm::InputTag src_; 

    edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
    edm::EDGetTokenT<reco::VertexCollection> offlinePrimaryVerticesToken_;
    edm::EDGetTokenT<SiPixelClusterShapeCache> pixelClusterShapeCacheToken_;
    edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> pixelSrcToken_;
    edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> tPixelDigi;
    //edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeometryToken_;
  };

  
  SiPixelPhase1TrackClusters::SiPixelPhase1TrackClusters(const edm::ParameterSet& iConfig)
      : SiPixelPhase1Base(iConfig), applyVertexCut_(iConfig.getUntrackedParameter<bool>("VertexCut", true)) {
    tracksToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));

    offlinePrimaryVerticesToken_ =
        applyVertexCut_ ? consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))
                        : edm::EDGetTokenT<reco::VertexCollection>();

    pixelClusterShapeCacheToken_ =
        consumes<SiPixelClusterShapeCache>(iConfig.getParameter<edm::InputTag>("clusterShapeCache"));

    pixelSrcToken_ = consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("clusters"));

    //src_ =  iConfighttps://github.com/cms-analysis/DPGAnalysis-SiPixelTools/blob/master/HitAnalyzer/test/PixDigisTest.cc#L1132.getParameter<edm::InputTag>( "src" );
    tPixelDigi =
      consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("src"));
  }

  void SiPixelPhase1TrackClusters::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    if (!checktrigger(iEvent, iSetup, DCS))
      return;

    if (histo.size() != ENUM_SIZE) {
      edm::LogError("SiPixelPhase1TrackClusters")
          << "incompatible configuration " << histo.size() << "!=" << ENUM_SIZE << std::endl;
      return;
    }

    // get geometry
    edm::ESHandle<TrackerGeometry> tracker;
    iSetup.get<TrackerDigiGeometryRecord>().get(tracker);
    assert(tracker.isValid());
    auto const& tracker_off = *tracker;

    edm::ESHandle<TrackerTopology> tTopoHandle;
    iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
    auto const& tkTpl = *tTopoHandle;

    edm::ESHandle<ClusterShapeHitFilter> shapeFilterH;
    iSetup.get<CkfComponentsRecord>().get("ClusterShapeHitFilter", shapeFilterH);
    auto const& shapeFilter = *shapeFilterH;

    edm::Handle<reco::VertexCollection> vertices;
    if (applyVertexCut_) {
      iEvent.getByToken(offlinePrimaryVerticesToken_, vertices);
      if (!vertices.isValid() || vertices->empty())
        return;
    }

    //get the map
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(tracksToken_, tracks);

    if (!tracks.isValid()) {
      edm::LogWarning("SiPixelPhase1TrackClusters") << "track collection is not valid";
      return;
    }

    edm::Handle<SiPixelClusterShapeCache> pixelClusterShapeCacheH;
    iEvent.getByToken(pixelClusterShapeCacheToken_, pixelClusterShapeCacheH);
    if (!pixelClusterShapeCacheH.isValid()) {
      edm::LogWarning("SiPixelPhase1TrackClusters") << "PixelClusterShapeCache collection is not valid";
      return;
    }
    auto const& pixelClusterShapeCache = *pixelClusterShapeCacheH;

    std::unordered_set<const SiPixelCluster*> rHSiPixelClusters;

    for (auto const& track : *tracks) {
      if (applyVertexCut_ &&
          (track.pt() < 0.75 || std::abs(track.dxy((*vertices)[0].position())) > 5 * track.dxyError()))
        continue;

      bool isBpixtrack = false, isFpixtrack = false, crossesPixVol = false;

      // find out whether track crosses pixel fiducial volume (for cosmic tracks)
      auto d0 = track.d0(), dz = track.dz();
      if (std::abs(d0) < 16 && std::abs(dz) < 50)
        crossesPixVol = true;

      auto etatk = track.eta();

      auto const& trajParams = track.extra()->trajParams();
      assert(trajParams.size() == track.recHitsSize());
      auto hb = track.recHitsBegin();

      for (unsigned int h = 0; h < track.recHitsSize(); h++) {
        auto hit = *(hb + h);
        if (!hit->isValid())
          continue;
        auto id = hit->geographicalId();

        // check that we are in the pixel
        auto subdetid = (id.subdetId());
        if (subdetid == PixelSubdetector::PixelBarrel)
          isBpixtrack = true;
        if (subdetid == PixelSubdetector::PixelEndcap)
          isFpixtrack = true;
        if (subdetid != PixelSubdetector::PixelBarrel && subdetid != PixelSubdetector::PixelEndcap)
          continue;
        bool iAmBarrel = subdetid == PixelSubdetector::PixelBarrel;

        // PXB_L4 IS IN THE OTHER WAY
        // CAN BE XORed BUT LETS KEEP THINGS SIMPLE
        bool iAmOuter = ((tkTpl.pxbLadder(id) % 2 == 1) && tkTpl.pxbLayer(id) != 4) ||
                        ((tkTpl.pxbLadder(id) % 2 != 1) && tkTpl.pxbLayer(id) == 4);

        auto pixhit = dynamic_cast<const SiPixelRecHit*>(hit->hit());
        if (!pixhit)
          continue;

        auto geomdetunit = dynamic_cast<const PixelGeomDetUnit*>(pixhit->detUnit());
        auto const& topol = geomdetunit->specificTopology();

        // get the cluster
        auto clustp = pixhit->cluster();
        if (clustp.isNull())
          continue;
        auto const& cluster = *clustp;

	const SiPixelCluster* rHcluster = &*(pixhit->cluster());

	rHSiPixelClusters.insert(rHcluster);
	/*
	edm::Handle< edm::DetSetVector<PixelDigi> > pixelDigis;
	iEvent.getByToken(tPixelDigi, pixelDigis);

	edm::DetSetVector<PixelDigi>::const_iterator DSViter;

	float total_digADC = 0.;
	*/
	//std::cout << "rHcluster inserted" << endl;
        const std::vector<SiPixelCluster::Pixel> pixelsVec = cluster.pixels();
        for (unsigned int i = 0; i < pixelsVec.size(); ++i) {
          float pixx = pixelsVec[i].x;  // index as float=iteger, row index
          float pixy = pixelsVec[i].y;  // same, col index

          bool bigInX = topol.isItBigPixelInX(int(pixx)); // dip solo da pixel
          bool bigInY = topol.isItBigPixelInY(int(pixy));// dip solo da pixel
          float pixel_charge = pixelsVec[i].adc;
	  /*
 	  for(DSViter = pixelDigis->begin(); DSViter != pixelDigis->end(); DSViter++) {
	    edm::DetSet<PixelDigi>::const_iterator di;
	    for(di = DSViter->data.begin(); di != DSViter->data.end(); di++) {
	      float digADC = di->adc(); // charge, modifued to unsiged short
	      float digx = di->row(); // row
	      float digy = di->column(); // column
	      if(digx == pixx && digy == pixy){
		total_digADC = total_digADC + digADC;
		cout << endl << "pix adc: " << pixel_charge << "\tdig adc: " << digADC << endl;
		//cout << "pix x: " << pixx << "\tdig x: " << digx << endl;
		//cout << "pix y: " << pixy << "\tdig y: " << digy << endl;
	      }
	      else continue;
	    }
	  }
	  cout << endl << "total digis adc so far: " << total_digADC << endl;
	  */

          if (bigInX == true || bigInY == true) {//si potrebbe fare a parte
            histo[ON_TRACK_BIGPIXELCHARGE].fill(pixel_charge, id, &iEvent); 
          } else {
            histo[ON_TRACK_NOTBIGPIXELCHARGE].fill(pixel_charge, id, &iEvent);
          }
        }  // End loop over pixels

	//float digiscluster_ratio = cluster.charge() / total_digADC;

	//histo[DIGIS_OVER_CLUSTER_TOTCHARGE].fill(digiscluster_ratio, id, &iEvent); 
	//histo[DIGIS_OVER_CLUSTER_TOTCHARGE_2D].fill(cluster.charge(), total_digADC, id, &iEvent); 

	//cout << endl << "cluster charge: " << cluster.charge() << "\ttotal digi adc: " << total_digADC << endl;

        auto const& ltp = trajParams[h];

        auto localDir = ltp.momentum() / ltp.momentum().mag();

        // correct charge for track impact angle
        auto charge = cluster.charge() * ltp.absdz();//per qualsiasi cluster
	
        auto clustgp = pixhit->globalPosition();  // from rechit

        int part;
        ClusterData::ArrayType meas;
        std::pair<float, float> pred;
        if (shapeFilter.getSizes(*pixhit, localDir, pixelClusterShapeCache, part, meas, pred)) {
          auto shape = shapeFilter.isCompatible(*pixhit, localDir, pixelClusterShapeCache);
          unsigned shapeVal = (shape ? 1 : 0);

          if (iAmBarrel) {
            if (iAmOuter) {
              histo[ON_TRACK_SIZE_X_OUTER].fill(pred.first, cluster.sizeX(), id, &iEvent);
              histo[ON_TRACK_SIZE_Y_OUTER].fill(pred.second, cluster.sizeY(), id, &iEvent);
              histo[ON_TRACK_SIZE_XY_OUTER].fill(cluster.sizeY(), cluster.sizeX(), id, &iEvent);

              histo[ON_TRACK_SHAPE_OUTER].fill(shapeVal, id, &iEvent);
            } else {
              histo[ON_TRACK_SIZE_X_INNER].fill(pred.first, cluster.sizeX(), id, &iEvent);
              histo[ON_TRACK_SIZE_Y_INNER].fill(pred.second, cluster.sizeY(), id, &iEvent);
              histo[ON_TRACK_SIZE_XY_INNER].fill(cluster.sizeY(), cluster.sizeX(), id, &iEvent);

              histo[ON_TRACK_SHAPE_INNER].fill(shapeVal, id, &iEvent);
            }
          } else {
            histo[ON_TRACK_SIZE_X_F].fill(pred.first, cluster.sizeX(), id, &iEvent);
            histo[ON_TRACK_SIZE_Y_F].fill(pred.second, cluster.sizeY(), id, &iEvent);
            histo[ON_TRACK_SIZE_XY_F].fill(cluster.sizeY(), cluster.sizeX(), id, &iEvent);
          }
          histo[ON_TRACK_SHAPE].fill(shapeVal, id, &iEvent);
        }

        for (int i = 0; i < cluster.size(); i++) {
          SiPixelCluster::Pixel const& vecipxl = cluster.pixel(i);
          histo[DIGIS_HITMAP_ON_TRACK].fill(id, &iEvent, vecipxl.y, vecipxl.x);
          histo[ON_TRACK_NDIGIS].fill(id, &iEvent);
        }

        histo[ON_TRACK_NCLUSTERS].fill(id, &iEvent);//qualsiasi cluster
        histo[ON_TRACK_CHARGE].fill(charge, id, &iEvent);//qualsiasi cluster
        histo[ON_TRACK_SIZE].fill(cluster.size(), id, &iEvent);//qualsiasi cluster

        histo[ON_TRACK_POSITIONB].fill(clustgp.z(), clustgp.phi(), id, &iEvent);//qualsiasi cluster
        histo[ON_TRACK_POSITIONF].fill(clustgp.x(), clustgp.y(), id, &iEvent);//qualsiasi cluster

        histo[CHARGE_VS_SIZE_ON_TRACK].fill(cluster.size(), charge, id, &iEvent);//qualsiasi cluster

        if (iAmBarrel)  // Avoid mistakes even if specification < should > handle it
        {
          if (iAmOuter) {
            histo[SIZE_VS_ETA_ON_TRACK_OUTER].fill(etatk, cluster.sizeY(), id, &iEvent);
            histo[ON_TRACK_CHARGE_OUTER].fill(charge, id, &iEvent);
          } else {
            histo[SIZE_VS_ETA_ON_TRACK_INNER].fill(etatk, cluster.sizeY(), id, &iEvent);
            histo[ON_TRACK_CHARGE_INNER].fill(charge, id, &iEvent);
          }
        }
      }

      // statistics on tracks
      histo[NTRACKS].fill(1, DetId(0), &iEvent);
      if (isBpixtrack || isFpixtrack)
        histo[NTRACKS].fill(2, DetId(0), &iEvent);
      if (isBpixtrack)
        histo[NTRACKS].fill(3, DetId(0), &iEvent);
      if (isFpixtrack)
        histo[NTRACKS].fill(4, DetId(0), &iEvent);

      if (crossesPixVol) {
        if (isBpixtrack || isFpixtrack)
          histo[NTRACKS_INVOLUME].fill(1, DetId(0), &iEvent);
        else
          histo[NTRACKS_INVOLUME].fill(0, DetId(0), &iEvent);
      }
    }

    
    //statistics clusters off tracks
    edm::Handle<edmNew::DetSetVector<SiPixelCluster>> inputPixel;
    iEvent.getByToken(pixelSrcToken_, inputPixel);

    if (!inputPixel.isValid())
      return;

    bool hasClusters = false;

    edmNew::DetSetVector<SiPixelCluster>::const_iterator it;    

    for (it = inputPixel->begin(); it != inputPixel->end(); ++it) {
      auto id = DetId(it->detId());

      const PixelGeomDetUnit* theGeomDet = dynamic_cast<const PixelGeomDetUnit*>(tracker_off.idToDet(id));
      const PixelTopology& topol = theGeomDet->specificTopology();

      for (SiPixelCluster const& cluster : *it) {
	if (!(rHSiPixelClusters.find(&cluster) == rHSiPixelClusters.end())) continue;
	std::cout << "Off track" << std::endl;

	int row = cluster.x() - 0.5, col = cluster.y() - 0.5;
	const std::vector<SiPixelCluster::Pixel> pixelsVec = cluster.pixels();
        for (unsigned int i = 0; i < pixelsVec.size(); ++i) {
          float pixx = pixelsVec[i].x;  // index as float=iteger, row index
          float pixy = pixelsVec[i].y;  // same, col index

	  cout << endl << i;
	  /*
	  if(i > 0){
	    cout << endl << "exit" << endl;
	    break;
	  }
	  */
          bool bigInX = topol.isItBigPixelInX(int(pixx)); // dip solo da pixel
          bool bigInY = topol.isItBigPixelInY(int(pixy));// dip solo da pixel
          float pixel_charge = pixelsVec[i].adc;

          if (bigInX == true || bigInY == true) {//si potrebbe fare a parte
            histo[OFF_TRACK_BIGPIXELCHARGE].fill(pixel_charge, id, &iEvent); 
	    std::cout << "BigPixel charge plotted" << std::endl;
          } 
	  else {
            histo[OFF_TRACK_NOTBIGPIXELCHARGE].fill(pixel_charge, id, &iEvent);
	    std::cout << "NotBigPixel charge plotted" << std::endl;
          }
        }  // End loop over pixels
	
	for (int i = 0; i < cluster.size(); i++) {
          SiPixelCluster::Pixel const& vecipxl = cluster.pixel(i);
          histo[DIGIS_HITMAP_OFF_TRACK].fill(id, &iEvent, vecipxl.y, vecipxl.x);
          histo[OFF_TRACK_NDIGIS].fill(id, &iEvent);
        }

	//histo[OFF_TRACK_READOUT_CHARGE].fill(double(cluster.charge()), id, &iEvent, col, row);
        histo[OFF_TRACK_CHARGE].fill(double(cluster.charge()), id, &iEvent, col, row);
        histo[OFF_TRACK_SIZE].fill(double(cluster.size()), id, &iEvent, col, row);
        histo[OFF_TRACK_SIZEX].fill(double(cluster.sizeX()), id, &iEvent, col, row);
        histo[OFF_TRACK_SIZEY].fill(double(cluster.sizeY()), id, &iEvent, col, row);
        histo[OFF_TRACK_NCLUSTERS].fill(id, &iEvent, col, row);

	/*
        if (cluster.size() > 1) {
          histo[OFF_TRACK_READOUT_NCLUSTERS].fill(id, &iEvent);
        }
	*/

        LocalPoint clustlp = topol.localPosition(MeasurementPoint(cluster.x(), cluster.y()));
        GlobalPoint clustgp = theGeomDet->surface().toGlobal(clustlp);
        
	histo[OFF_TRACK_POSITION_B].fill(clustgp.z(), clustgp.phi(), id, &iEvent);
        histo[OFF_TRACK_POSITION_F].fill(clustgp.x(), clustgp.y(), id, &iEvent);
        histo[OFF_TRACK_POSITION_XZ].fill(clustgp.x(), clustgp.z(), id, &iEvent);
        histo[OFF_TRACK_POSITION_YZ].fill(clustgp.y(), clustgp.z(), id, &iEvent);
        histo[SIZE_VS_ETA_OFF_TRACK].fill(clustgp.eta(), cluster.sizeY(), id, &iEvent);
	histo[CHARGE_VS_ETA_OFF_TRACK].fill(clustgp.eta(), cluster.charge(), id, &iEvent);
	histo[CHARGE_VS_SIZE_OFF_TRACK].fill(cluster.sizeY(), cluster.charge(), id, &iEvent);

      }

    }
    
    histo[ON_TRACK_NCLUSTERS].executePerEventHarvesting(&iEvent);
    histo[ON_TRACK_NDIGIS].executePerEventHarvesting(&iEvent);
    histo[OFF_TRACK_NCLUSTERS].executePerEventHarvesting(&iEvent);
    histo[OFF_TRACK_NDIGIS].executePerEventHarvesting(&iEvent);
  }

}  // namespace

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SiPixelPhase1TrackClusters);
