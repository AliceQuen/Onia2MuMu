// ========== 阶段1: Jpsi候选预缓存 ==========
  // 通过muon双重循环构建所有Jpsi候选并缓存，避免重复拟合
  std::vector<JpsiCandidate> jpsiCandidates;
  jpsiCandidates.reserve(muPlus.size() * muMinus.size());
  
  // 记录每个muon的触发匹配状态
  std::map<edm::View<pat::Muon>::const_iterator, bool> muonFilterMatchMap;
  
  for (const auto& muPlusIter : muPlus) {
    TrackRef muTrack1 = muPlusIter->track();
    
    for (const auto& muMinusIter : muMinus) {
      TrackRef muTrack2 = muMinusIter->track();
      
      // 质量预筛选（1-4.5 GeV）
      double invMass = (muPlusIter->p4() + muMinusIter->p4()).mass();
      if (!(1.0 < invMass && invMass < 4.5)) {
        #if DEBUG == 2
        continue_counts_[__LINE__]++;
        #endif
        continue;
      }
      
      TransientTrack muonTT1(muTrack1, &(bFieldHandle));
      TransientTrack muonTT2(muTrack2, &(bFieldHandle));
      
      KinematicParticleFactoryFromTransientTrack pmumuFactory;
      ParticleMass muon_mass = MU_MASS;
      float muon_sigma = MU_MASSERR;
      float chi = 0.;
      float ndf = 0.;
      
      vector<RefCountedKinematicParticle> muonParticles;
      muonParticles.push_back(pmumuFactory.particle(muonTT1, muon_mass, chi, ndf, muon_sigma));
      muonParticles.push_back(pmumuFactory.particle(muonTT2, muon_mass, chi, ndf, muon_sigma));
      
      KinematicParticleVertexFitter fitter;
      RefCountedKinematicTree vertexFitTree;
      Error_t = false;
      try {
        vertexFitTree = fitter.fit(muonParticles);
      } catch (...) {
        Error_t = true;
      }
      
      if (Error_t || !vertexFitTree->isValid()) {
        #if DEBUG == 2
        continue_counts_[__LINE__]++;
        total_continue_++;
        #endif
        continue;
      }
      
      vertexFitTree->movePointerToTheTop();
      RefCountedKinematicParticle jpsiParticle = vertexFitTree->currentParticle();
      RefCountedKinematicVertex jpsiVertex = vertexFitTree->currentDecayVertex();
      
      double vtxProb = ChiSquaredProbability(
          (double)(jpsiVertex->chiSquared()),
          (double)(jpsiVertex->degreesOfFreedom()));
      
      // 顶点概率筛选（>1%）
      if (vtxProb < jpsi_vtxprob_cut) {
        #if DEBUG == 2
        continue_counts_[__LINE__]++;
        total_continue_++;
        #endif
        continue;
      }
      
      // 质量窗口筛选（Jpsi: ±150 MeV）
      double jpsiMass = jpsiParticle->currentState().mass();
      if (fabs(jpsiMass - jpsi_nominal_mass) > jpsi_mass_window) {
        #if DEBUG == 2
        continue_counts_[__LINE__]++;
        total_continue_++;
        #endif
        continue;
      }
      
      // 获取触发匹配状态
      bool matchPlus = false, matchMinus = false;
      for (unsigned int JpsiFilt = 0; JpsiFilt < FiltersForJpsi_.size(); JpsiFilt++) {
        if (hltresults.isValid()) {
          for (auto i = muPlusIter->triggerObjectMatches().begin();
               i != muPlusIter->triggerObjectMatches().end(); i++) {
            pat::TriggerObjectStandAlone tempTriggerObject(*i);
            tempTriggerObject.unpackFilterLabels(iEvent, *hltresults);
            if (tempTriggerObject.hasFilterLabel(FiltersForJpsi_[JpsiFilt])) {
              matchPlus = true;
              break;
            }
          }
          if (matchPlus) break;
        }
      }
      
      for (unsigned int JpsiFilt = 0; JpsiFilt < FiltersForJpsi_.size(); JpsiFilt++) {
        if (hltresults.isValid()) {
          for (auto i = muMinusIter->triggerObjectMatches().begin();
               i != muMinusIter->triggerObjectMatches().end(); i++) {
            pat::TriggerObjectStandAlone tempTriggerObject(*i);
            tempTriggerObject.unpackFilterLabels(iEvent, *hltresults);
            if (tempTriggerObject.hasFilterLabel(FiltersForJpsi_[JpsiFilt])) {
              matchMinus = true;
              break;
            }
          }
          if (matchMinus) break;
        }
      }
      
      // 缓存Jpsi候选
      JpsiCandidate candidate;
      candidate.muPlus = muPlusIter;
      candidate.muMinus = muMinusIter;
      candidate.mass = jpsiMass;
      candidate.vtxProb = vtxProb;
      candidate.massErr = (jpsiParticle->currentState().kinematicParametersError().matrix()(6, 6) > 0) 
                         ? sqrt(jpsiParticle->currentState().kinematicParametersError().matrix()(6, 6)) 
                         : -9;
      candidate.p4.SetPtEtaPhiM(muPlusIter->track()->pt(), muPlusIter->track()->eta(),
                                muPlusIter->track()->phi(), MU_MASS);
      candidate.p4 += TLorentzVector(muMinusIter->track()->pt(), muMinusIter->track()->eta(),
                                     muMinusIter->track()->phi(), MU_MASS);
      candidate.kinematicParticle = jpsiParticle;
      candidate.vertex = jpsiVertex;
      candidate.muonTT1 = muonTT1;
      candidate.muonTT2 = muonTT2;
      candidate.filterMatchPlus = matchPlus;
      candidate.filterMatchMinus = matchMinus;
      
      jpsiCandidates.push_back(candidate);
    }
  }
  
  #if DEBUG == 1
  std::cout << "[DEBUG] Number of Jpsi candidates: " << jpsiCandidates.size() << std::endl;
  #endif
  
  // 检查是否有足够的Jpsi候选
  if (jpsiCandidates.size() < 2) {
    #if DEBUG == 1
    std::cout << "[DEBUG] Not enough Jpsi candidates, skipping event" << std::endl;
    #endif
    return;
  }

  // ========== 阶段2: Jpsi候选两两组合 ==========
  // 从缓存中选择两个不同的Jpsi进行组合
  for (size_t i = 0; i < jpsiCandidates.size(); ++i) {
    const JpsiCandidate& jpsi1 = jpsiCandidates[i];
    
    for (size_t j = i + 1; j < jpsiCandidates.size(); ++j) {
      const JpsiCandidate& jpsi2 = jpsiCandidates[j];
      
      // 检查是否使用了相同的muon
      if (jpsi1.muPlus == jpsi2.muPlus || jpsi1.muPlus == jpsi2.muMinus ||
          jpsi1.muMinus == jpsi2.muPlus || jpsi1.muMinus == jpsi2.muMinus) {
        #if DEBUG == 2
        continue_counts_[__LINE__]++;
        #endif
        continue;
      }
      
      // 触发匹配检查：jpsi1的两个muon匹配 或 jpsi2的两个muon匹配
      bool passTrigger = (jpsi1.filterMatchPlus && jpsi1.filterMatchMinus) ||
                         (jpsi2.filterMatchPlus && jpsi2.filterMatchMinus);
      if (!passTrigger) {
        #if DEBUG == 2
        continue_counts_[__LINE__]++;
        #endif
        continue;
      }

      // ========== 阶段3: 与Track组合 ==========
      for (const auto& trackPlusIter : trackPlus) {
        if (!trackPlusIter->hasTrackDetails() || trackPlusIter->charge() == 0) {
          #if DEBUG == 2
          continue_counts_[__LINE__]++;
          #endif
          continue;
        }
        const reco::Track* track1 = trackPlusIter->bestTrack();
        if (!track1) {
          #if DEBUG == 2
          continue_counts_[__LINE__]++;
          #endif
          continue;
        }

        for (const auto& trackMinusIter : trackMinus) {
          if (!trackMinusIter->hasTrackDetails() || trackMinusIter->charge() == 0) {
            #if DEBUG == 2
            continue_counts_[__LINE__]++;
            #endif
            continue;
          }
          const reco::Track* track2 = trackMinusIter->bestTrack();
          if (!track2) {
            #if DEBUG == 2
            continue_counts_[__LINE__]++;
            #endif
            continue;
          }

          TLorentzVector P4_Track1, P4_Track2, P4_Jpsipipi;
          P4_Track1.SetPtEtaPhiM(trackPlusIter->pt(), trackPlusIter->eta(),
                                 trackPlusIter->phi(), PI_MASS);
          P4_Track2.SetPtEtaPhiM(trackMinusIter->pt(), trackMinusIter->eta(),
                                 trackMinusIter->phi(), PI_MASS);
          P4_Jpsipipi = jpsi1.p4 + P4_Track1 + P4_Track2;

          if (P4_Track1.DeltaR(P4_Jpsipipi) > pionDRcut) {
            #if DEBUG == 2
            continue_counts_[__LINE__]++;
            total_continue_++;
            #endif
            continue;
          }
          if (P4_Track2.DeltaR(P4_Jpsipipi) > pionDRcut) {
            #if DEBUG == 2
            continue_counts_[__LINE__]++;
            total_continue_++;
            #endif
            continue;
          }

          TransientTrack trackTT1(*track1, &(bFieldHandle));
          TransientTrack trackTT2(*track2, &(bFieldHandle));
          KinematicParticleFactoryFromTransientTrack JPiPiFactory;
          ParticleMass pion_mass = PI_MASS;
          float pion_sigma = PI_MASSERR;
          float chi = 0.;
          float ndf = 0.;

          // ========== 四粒子带J/psi质量约束拟合 ==========
          vector<RefCountedKinematicParticle> JPiPiParticles;
          JPiPiParticles.push_back(JPiPiFactory.particle(
              trackTT1, pion_mass, chi, ndf, pion_sigma));
          JPiPiParticles.push_back(JPiPiFactory.particle(
              trackTT2, pion_mass, chi, ndf, pion_sigma));
          JPiPiParticles.push_back(JPiPiFactory.particle(
              jpsi1.muonTT1, MU_MASS, chi, ndf, MU_MASSERR));
          JPiPiParticles.push_back(JPiPiFactory.particle(
              jpsi1.muonTT2, MU_MASS, chi, ndf, MU_MASSERR));

          KinematicConstraint* jpsiMassConstraint = new MassKinematicConstraint(JPSI_MASS_NOMINAL, 2, 3);
          vector<KinematicConstraint*> JPiPiConstraints;
          JPiPiConstraints.push_back(jpsiMassConstraint);

          KinematicConstrainedVertexFitter JPiPi_fitter;
          RefCountedKinematicTree JPiPiVertexFitTree;
          Error_t = false;
          try {
            JPiPiVertexFitTree = JPiPi_fitter.fit(JPiPiParticles, JPiPiConstraints);
          } catch (...) {
            Error_t = true;
          }
          if (Error_t || !(JPiPiVertexFitTree->isValid())) {
            delete jpsiMassConstraint;
            continue;
          }
          
          JPiPiVertexFitTree->movePointerToTheTop();
          RefCountedKinematicParticle JPiPi_vFit_constrained =
              JPiPiVertexFitTree->currentParticle();
          RefCountedKinematicVertex JPiPi_vFit_vertex_constrained =
              JPiPiVertexFitTree->currentDecayVertex();

          double JPiPi_vtxprob = ChiSquaredProbability(
              (double)(JPiPi_vFit_vertex_constrained->chiSquared()),
              (double)(JPiPi_vFit_vertex_constrained->degreesOfFreedom()));
          
          if (JPiPi_vFit_constrained->currentState().mass() > 4.5) {
            delete jpsiMassConstraint;
            continue;
          }
          
          if (JPiPi_vtxprob < psi2s_vtxprob_cut) {
            delete jpsiMassConstraint;
            continue;
          }
          
          double px = JPiPi_vFit_constrained->currentState().kinematicParameters().momentum().x();
          double py = JPiPi_vFit_constrained->currentState().kinematicParameters().momentum().y();
          double psi2s_pt = sqrt(px*px + py*py);
          if (psi2s_pt <= psi2s_pt_cut) {
            delete jpsiMassConstraint;
            continue;
          }

          delete jpsiMassConstraint;

          // ========== 保存Psi2S拟合结果 ==========
          Psi2S_mass = JPiPi_vFit_constrained->currentState().mass();
          Psi2S_VtxProb = JPiPi_vtxprob;
          Psi2S_px = JPiPi_vFit_constrained->currentState()
                                     .kinematicParameters()
                                     .momentum()
                                     .x();
          Psi2S_py = JPiPi_vFit_constrained->currentState()
                                     .kinematicParameters()
                                     .momentum()
                                     .y();
          Psi2S_pz = JPiPi_vFit_constrained->currentState()
                                     .kinematicParameters()
                                     .momentum()
                                     .z();
          if (JPiPi_vFit_constrained->currentState()
                  .kinematicParametersError()
                  .matrix()(6, 6) > 0) {
            Psi2S_massErr = sqrt(JPiPi_vFit_constrained->currentState()
                                        .kinematicParametersError()
                                        .matrix()(6, 6));
          } else {
            Psi2S_massErr = -9;
          }

          ROOT::Math::PxPyPzMVector Psi2S_vec(Psi2S_px, Psi2S_py, Psi2S_pz, Psi2S_mass);
          Psi2S_pt = Psi2S_vec.Pt();
          Psi2S_absEta = fabs(Psi2S_vec.Eta());

          // ========== 六粒子三种假设质量约束拟合 ==========
          vector<RefCountedKinematicParticle> X_Particles;
          X_Particles.push_back(JPiPiFactory.particle(
              trackTT1, pion_mass, chi, ndf, pion_sigma));
          X_Particles.push_back(JPiPiFactory.particle(
              trackTT2, pion_mass, chi, ndf, pion_sigma));
          X_Particles.push_back(JPiPiFactory.particle(
              jpsi1.muonTT1, MU_MASS, chi, ndf, MU_MASSERR));
          X_Particles.push_back(JPiPiFactory.particle(
              jpsi1.muonTT2, MU_MASS, chi, ndf, MU_MASSERR));
          X_Particles.push_back(JPiPiFactory.particle(
              jpsi2.muonTT1, MU_MASS, chi, ndf, MU_MASSERR));
          X_Particles.push_back(JPiPiFactory.particle(
              jpsi2.muonTT2, MU_MASS, chi, ndf, MU_MASSERR));

          // ========== 假设PJ: mumupipi(Jpsi) + mumup(Jpsi) ==========
          {
            vector<KinematicConstraint*> PJConstraints;
            KinematicConstraint* JpsiConstraint1 = new MassKinematicConstraint(JPSI_MASS_NOMINAL, 2, 3);
            KinematicConstraint* JpsiConstraint2 = new MassKinematicConstraint(JPSI_MASS_NOMINAL, 4, 5);
            PJConstraints.push_back(JpsiConstraint1);
            PJConstraints.push_back(JpsiConstraint2);

            KinematicConstrainedVertexFitter PJ_fitter;
            RefCountedKinematicTree PJ_VertexFitTree;
            bool PJ_Error = false;
            try {
              PJ_VertexFitTree = PJ_fitter.fit(X_Particles, PJConstraints);
            } catch (...) {
              PJ_Error = true;
            }
            if (PJ_Error || !(PJ_VertexFitTree->isValid())) {
              delete JpsiConstraint1;
              delete JpsiConstraint2;
            } else {
              PJ_VertexFitTree->movePointerToTheTop();
              RefCountedKinematicParticle PJ_vFit = PJ_VertexFitTree->currentParticle();
              RefCountedKinematicVertex PJ_vFit_vertex = PJ_VertexFitTree->currentDecayVertex();
              KinematicParameters PJ_kPara = PJ_vFit->currentState().kinematicParameters();

              X_PJ_mass = PJ_vFit->currentState().mass();
              X_PJ_VtxProb = ChiSquaredProbability(
                  (double)(PJ_vFit_vertex->chiSquared()),
                  (double)(PJ_vFit_vertex->degreesOfFreedom()));
              X_PJ_px = PJ_kPara.momentum().x();
              X_PJ_py = PJ_kPara.momentum().y();
              X_PJ_pz = PJ_kPara.momentum().z();
              if (PJ_vFit->currentState().kinematicParametersError().matrix()(6, 6) > 0) {
                X_PJ_massErr = sqrt(PJ_vFit->currentState().kinematicParametersError().matrix()(6, 6));
              } else {
                X_PJ_massErr = -9;
              }
              ROOT::Math::PxPyPzMVector PJ_vec(X_PJ_px, X_PJ_py, X_PJ_pz, X_PJ_mass);
              X_PJ_pt = PJ_vec.Pt();
              X_PJ_absEta = fabs(PJ_vec.Eta());
              
              delete JpsiConstraint1;
              delete JpsiConstraint2;
            }
          }

          // ========== 假设XJ: mumupipi(X3872) + mumup(Jpsi) ==========
          {
            vector<KinematicConstraint*> XJConstraints;
            vector<int> X3872_particleList;
            X3872_particleList.push_back(0);
            X3872_particleList.push_back(1);
            X3872_particleList.push_back(2);
            X3872_particleList.push_back(3);
            KinematicConstraint* X3872Constraint = new MassKinematicConstraint(X3872_MASS_NOMINAL, X3872_particleList);
            KinematicConstraint* JpsiConstraint2 = new MassKinematicConstraint(JPSI_MASS_NOMINAL, 4, 5);
            XJConstraints.push_back(X3872Constraint);
            XJConstraints.push_back(JpsiConstraint2);

            KinematicConstrainedVertexFitter XJ_fitter;
            RefCountedKinematicTree XJ_VertexFitTree;
            bool XJ_Error = false;
            try {
              XJ_VertexFitTree = XJ_fitter.fit(X_Particles, XJConstraints);
            } catch (...) {
              XJ_Error = true;
            }
            if (XJ_Error || !(XJ_VertexFitTree->isValid())) {
              delete X3872Constraint;
              delete JpsiConstraint2;
            } else {
              XJ_VertexFitTree->movePointerToTheTop();
              RefCountedKinematicParticle XJ_vFit = XJ_VertexFitTree->currentParticle();
              RefCountedKinematicVertex XJ_vFit_vertex = XJ_VertexFitTree->currentDecayVertex();
              KinematicParameters XJ_kPara = XJ_vFit->currentState().kinematicParameters();

              X_XJ_mass = XJ_vFit->currentState().mass();
              X_XJ_VtxProb = ChiSquaredProbability(
                  (double)(XJ_vFit_vertex->chiSquared()),
                  (double)(XJ_vFit_vertex->degreesOfFreedom()));
              X_XJ_px = XJ_kPara.momentum().x();
              X_XJ_py = XJ_kPara.momentum().y();
              X_XJ_pz = XJ_kPara.momentum().z();
              if (XJ_vFit->currentState().kinematicParametersError().matrix()(6, 6) > 0) {
                X_XJ_massErr = sqrt(XJ_vFit->currentState().kinematicParametersError().matrix()(6, 6));
              } else {
                X_XJ_massErr = -9;
              }
              ROOT::Math::PxPyPzMVector XJ_vec(X_XJ_px, X_XJ_py, X_XJ_pz, X_XJ_mass);
              X_XJ_pt = XJ_vec.Pt();
              X_XJ_absEta = fabs(XJ_vec.Eta());
              
              delete X3872Constraint;
              delete JpsiConstraint2;
            }
          }

          // ========== 假设PP: mumupipi(Jpsi) + mumup(Psi2S) ==========
          {
            vector<KinematicConstraint*> PPConstraints;
            KinematicConstraint* JpsiConstraint1 = new MassKinematicConstraint(JPSI_MASS_NOMINAL, 2, 3);
            KinematicConstraint* Psi2SConstraint2 = new MassKinematicConstraint(PSI2S_MASS_NOMINAL, 4, 5);
            PPConstraints.push_back(JpsiConstraint1);
            PPConstraints.push_back(Psi2SConstraint2);

            KinematicConstrainedVertexFitter PP_fitter;
            RefCountedKinematicTree PP_VertexFitTree;
            bool PP_Error = false;
            try {
              PP_VertexFitTree = PP_fitter.fit(X_Particles, PPConstraints);
            } catch (...) {
              PP_Error = true;
            }
            if (PP_Error || !(PP_VertexFitTree->isValid())) {
              delete JpsiConstraint1;
              delete Psi2SConstraint2;
            } else {
              PP_VertexFitTree->movePointerToTheTop();
              RefCountedKinematicParticle PP_vFit = PP_VertexFitTree->currentParticle();
              RefCountedKinematicVertex PP_vFit_vertex = PP_VertexFitTree->currentDecayVertex();
              KinematicParameters PP_kPara = PP_vFit->currentState().kinematicParameters();

              X_PP_mass = PP_vFit->currentState().mass();
              X_PP_VtxProb = ChiSquaredProbability(
                  (double)(PP_vFit_vertex->chiSquared()),
                  (double)(PP_vFit_vertex->degreesOfFreedom()));
              X_PP_px = PP_kPara.momentum().x();
              X_PP_py = PP_kPara.momentum().y();
              X_PP_pz = PP_kPara.momentum().z();
              if (PP_vFit->currentState().kinematicParametersError().matrix()(6, 6) > 0) {
                X_PP_massErr = sqrt(PP_vFit->currentState().kinematicParametersError().matrix()(6, 6));
              } else {
                X_PP_massErr = -9;
              }
              ROOT::Math::PxPyPzMVector PP_vec(X_PP_px, X_PP_py, X_PP_pz, X_PP_mass);
              X_PP_pt = PP_vec.Pt();
              X_PP_absEta = fabs(PP_vec.Eta());
              
              delete JpsiConstraint1;
              delete Psi2SConstraint2;
            }
          }

          // ========== 保存Jpsi变量 ==========
          Jpsi1_mass = jpsi1.mass;
          Jpsi1_VtxProb = jpsi1.vtxProb;
          Jpsi1_massErr = jpsi1.massErr;
          ROOT::Math::PxPyPzMVector Jpsi1_vec(jpsi1.p4.Px(), jpsi1.p4.Py(), jpsi1.p4.Pz(), jpsi1.mass);
          Jpsi1_pt = Jpsi1_vec.Pt();
          Jpsi1_absEta = fabs(Jpsi1_vec.Eta());

          Jpsi2_mass = jpsi2.mass;
          Jpsi2_VtxProb = jpsi2.vtxProb;
          Jpsi2_massErr = jpsi2.massErr;
          ROOT::Math::PxPyPzMVector Jpsi2_vec(jpsi2.p4.Px(), jpsi2.p4.Py(), jpsi2.p4.Pz(), jpsi2.mass);
          Jpsi2_pt = Jpsi2_vec.Pt();
          Jpsi2_absEta = fabs(Jpsi2_vec.Eta());

          // ========== 保存muon变量 ==========
          const auto& iMuon1 = *jpsi1.muPlus;
          const auto& iMuon2 = *jpsi1.muMinus;
          const auto& iMuon3 = *jpsi2.muPlus;
          const auto& iMuon4 = *jpsi2.muMinus;

          mu1_hasFilterMatch = jpsi1.filterMatchPlus ? 1 : 0;
          mu2_hasFilterMatch = jpsi1.filterMatchMinus ? 1 : 0;
          mu3_hasFilterMatch = jpsi2.filterMatchPlus ? 1 : 0;
          mu4_hasFilterMatch = jpsi2.filterMatchMinus ? 1 : 0;

          mu1_px = iMuon1.px();
          mu1_py = iMuon1.py();
          mu1_pz = iMuon1.pz();
          mu1_pt = iMuon1.pt();
          mu1_absEta = fabs(iMuon1.eta());
          mu1_trackIso = iMuon1.trackIso();
          mu1_d0BS = iMuon1.dB(pat::Muon::BS2D);
          mu1_d0EBS = iMuon1.edB(pat::Muon::BS2D);
          mu1_d3dBS = iMuon1.dB(pat::Muon::BS3D);
          mu1_d3dEBS = iMuon1.edB(pat::Muon::BS3D);
          mu1_d0PV = iMuon1.dB(pat::Muon::PV2D);
          mu1_d0EPV = iMuon1.edB(pat::Muon::PV2D);
          mu1_dzPV = iMuon1.dB(pat::Muon::PVDZ);
          mu1_dzEPV = iMuon1.edB(pat::Muon::PVDZ);
          mu1_charge = iMuon1.charge();

          mu2_px = iMuon2.px();
          mu2_py = iMuon2.py();
          mu2_pz = iMuon2.pz();
          mu2_pt = iMuon2.pt();
          mu2_absEta = fabs(iMuon2.eta());
          mu2_trackIso = iMuon2.trackIso();
          mu2_d0BS = iMuon2.dB(pat::Muon::BS2D);
          mu2_d0EBS = iMuon2.edB(pat::Muon::BS2D);
          mu2_d3dBS = iMuon2.dB(pat::Muon::BS3D);
          mu2_d3dEBS = iMuon2.edB(pat::Muon::BS3D);
          mu2_d0PV = iMuon2.dB(pat::Muon::PV2D);
          mu2_d0EPV = iMuon2.edB(pat::Muon::PV2D);
          mu2_dzPV = iMuon2.dB(pat::Muon::PVDZ);
          mu2_dzEPV = iMuon2.edB(pat::Muon::PVDZ);
          mu2_charge = iMuon2.charge();

          mu3_px = iMuon3.px();
          mu3_py = iMuon3.py();
          mu3_pz = iMuon3.pz();
          mu3_pt = iMuon3.pt();
          mu3_absEta = fabs(iMuon3.eta());
          mu3_trackIso = iMuon3.trackIso();
          mu3_d0BS = iMuon3.dB(pat::Muon::BS2D);
          mu3_d0EBS = iMuon3.edB(pat::Muon::BS2D);
          mu3_d3dBS = iMuon3.dB(pat::Muon::BS3D);
          mu3_d3dEBS = iMuon3.edB(pat::Muon::BS3D);
          mu3_d0PV = iMuon3.dB(pat::Muon::PV2D);
          mu3_d0EPV = iMuon3.edB(pat::Muon::PV2D);
          mu3_dzPV = iMuon3.dB(pat::Muon::PVDZ);
          mu3_dzEPV = iMuon3.edB(pat::Muon::PVDZ);
          mu3_charge = iMuon3.charge();

          mu4_px = iMuon4.px();
          mu4_py = iMuon4.py();
          mu4_pz = iMuon4.pz();
          mu4_pt = iMuon4.pt();
          mu4_absEta = fabs(iMuon4.eta());
          mu4_trackIso = iMuon4.trackIso();
          mu4_d0BS = iMuon4.dB(pat::Muon::BS2D);
          mu4_d0EBS = iMuon4.edB(pat::Muon::BS2D);
          mu4_d3dBS = iMuon4.dB(pat::Muon::BS3D);
          mu4_d3dEBS = iMuon4.edB(pat::Muon::BS3D);
          mu4_d0PV = iMuon4.dB(pat::Muon::PV2D);
          mu4_d0EPV = iMuon4.edB(pat::Muon::PV2D);
          mu4_dzPV = iMuon4.dB(pat::Muon::PVDZ);
          mu4_dzEPV = iMuon4.edB(pat::Muon::PVDZ);
          mu4_charge = iMuon4.charge();

          // Count muon IDs for all 4 muons
          nLooseMuons = 0;
          nTightMuons = 0;
          nSoftMuons = 0;
          nMediumMuons = 0;

          if (iMuon1.isLooseMuon()) nLooseMuons++;
          if (iMuon1.isTightMuon(thePrimaryV)) nTightMuons++;
          if (iMuon1.isSoftMuon(thePrimaryV)) nSoftMuons++;
          if (iMuon1.isMediumMuon()) nMediumMuons++;

          if (iMuon2.isLooseMuon()) nLooseMuons++;
          if (iMuon2.isTightMuon(thePrimaryV)) nTightMuons++;
          if (iMuon2.isSoftMuon(thePrimaryV)) nSoftMuons++;
          if (iMuon2.isMediumMuon()) nMediumMuons++;

          if (iMuon3.isLooseMuon()) nLooseMuons++;
          if (iMuon3.isTightMuon(thePrimaryV)) nTightMuons++;
          if (iMuon3.isSoftMuon(thePrimaryV)) nSoftMuons++;
          if (iMuon3.isMediumMuon()) nMediumMuons++;

          if (iMuon4.isLooseMuon()) nLooseMuons++;
          if (iMuon4.isTightMuon(thePrimaryV)) nTightMuons++;
          if (iMuon4.isSoftMuon(thePrimaryV)) nSoftMuons++;
          if (iMuon4.isMediumMuon()) nMediumMuons++;

          // ========== 保存pion变量 ==========
          pi1_px = trackTT1.impactPointState().globalMomentum().x();
          pi1_py = trackTT1.impactPointState().globalMomentum().y();
          pi1_pz = trackTT1.impactPointState().globalMomentum().z();
          ROOT::Math::PxPyPzMVector pi1_vec(pi1_px, pi1_py, pi1_pz, PI_MASS);
          pi1_pt = pi1_vec.Pt();
          pi1_absEta = fabs(pi1_vec.Eta());

          pi2_px = trackTT2.impactPointState().globalMomentum().x();
          pi2_py = trackTT2.impactPointState().globalMomentum().y();
          pi2_pz = trackTT2.impactPointState().globalMomentum().z();
          ROOT::Math::PxPyPzMVector pi2_vec(pi2_px, pi2_py, pi2_pz, PI_MASS);
          pi2_pt = pi2_vec.Pt();
          pi2_absEta = fabs(pi2_vec.Eta());

          // ========== 计算DeltaR变量 ==========
          ROOT::Math::PxPyPzMVector mu1_vec(mu1_px, mu1_py, mu1_pz, MU_MASS);
          ROOT::Math::PxPyPzMVector mu2_vec(mu2_px, mu2_py, mu2_pz, MU_MASS);
          ROOT::Math::PxPyPzMVector mu3_vec(mu3_px, mu3_py, mu3_pz, MU_MASS);
          ROOT::Math::PxPyPzMVector mu4_vec(mu4_px, mu4_py, mu4_pz, MU_MASS);

          dR_mu1_mu2 = ROOT::Math::VectorUtil::DeltaR(mu1_vec, mu2_vec);
          dR_mu3_mu4 = ROOT::Math::VectorUtil::DeltaR(mu3_vec, mu4_vec);
          dR_pi1_pi2 = ROOT::Math::VectorUtil::DeltaR(pi1_vec, pi2_vec);
          dR_Psi2S_Jpsi1 = ROOT::Math::VectorUtil::DeltaR(Psi2S_vec, Jpsi1_vec);
          dR_Psi2S_Jpsi2 = ROOT::Math::VectorUtil::DeltaR(Psi2S_vec, Jpsi2_vec);
          dR_Psi2S_pi1 = ROOT::Math::VectorUtil::DeltaR(Psi2S_vec, pi1_vec);
          dR_Psi2S_pi2 = ROOT::Math::VectorUtil::DeltaR(Psi2S_vec, pi2_vec);

          X_One_Tree_->Fill();
        } // end trackMinus loop
      } // end trackPlus loop
    } // end jpsi2 loop
  } // end jpsi1 loop