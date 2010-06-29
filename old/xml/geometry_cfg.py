import FWCore.ParameterSet.Config as cms

from FWCore.MessageLogger.MessageLogger_cfi import *
process = cms.Process("IGUANA")
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.VisConfigurationService = cms.Service("VisConfigurationService",
					      Views = cms.untracked.vstring('3D Window'),
					      ContentProxies = cms.untracked.vstring('Simulation/Core',
										     'Simulation/Geometry',
										     'Reco/CMS Magnetic Field')
					     )

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

process.prod = cms.EDProducer("GeometryProducer",
			      UseMagneticField = cms.bool(False),
			      UseSensitiveDetectors = cms.bool(False),
			      MagneticField = cms.PSet(delta = cms.double(1.0))
			     )

process.p = cms.Path(process.prod)
