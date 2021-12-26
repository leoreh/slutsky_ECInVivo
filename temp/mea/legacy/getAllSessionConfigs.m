
function sessionConfigs = getAllSessionConfigs()
    sessionConfigs = struct;
   % sessionConfigs.bac1a = newConfig('1uM-bac/2.9.12_1uM_bac.mat',3);
    sessionConfigs.bac1b = newConfig('1uM-bac/13.3.13_1uM_bac.mat',4);
    sessionConfigs.bac1c = newConfig('1uM-bac/6.12.12_1uM_bac.mat',3);
    sessionConfigs.bac1d = newConfig('1uM-bac/19.11.12_1uM_bac.mat',3);
    sessionConfigs.bac1e = newConfig('1uM-bac/29.10.13_1uM_bac.mat',4);
    
    sessionConfigs.bac10a = newConfig('10uM-bac/16.3.13_10uM_bac.mat',4);
    sessionConfigs.bac10b = newConfig('10uM-bac/16.10.12_10uM_bac.mat',3);
    sessionConfigs.bac10c = newConfig('10uM -bac/4.12.12_10uM_bac.mat',3);
    sessionConfigs.bac10d = newConfig('10uM-bac/9.10.12_10uM_bac.mat',3);
    sessionConfigs.bac10e = newConfig('10uM-bac/3.10.13_10uM_bac.mat',4); % done
    sessionConfigs.bac10f = newConfig('10uM-bac/23.10.13_10uM_bac.mat',4);
    sessionConfigs.bac10g = newConfig('10uM-bac/23.11.13_10uM_bac.mat',4);
    
    sessionConfigs.cgp1a = newConfig('1uM-CGP/12.11.13_1uM_cgp.mat',4); % done
    sessionConfigs.cgp1b = newConfig('1uM-CGP/20.10.13_1uM_cgp.mat',3); % done
    sessionConfigs.cgp1c = newConfig('1uM-CGP/27.10.13_1uM_cgp.mat',4); % done
    
   % sessionConfigs.cnqx40a = newConfig('40uM-CNQX/10.2.13_40uM_CNQX.mat',4);
    sessionConfigs.cnqx40b = newConfig('40uM-CNQX/27.1.13_40uM_CNQX.mat',3);
    sessionConfigs.cnqx40c = newConfig('40uM-CNQX/27.3.13_40_uM_CNQX.mat',4);
    sessionConfigs.cnqx40d = newConfig('40uM-CNQX/9.5.13_10uM_CNQX.mat',4); % done
    sessionConfigs.cnqx40e = newConfig('40uM-CNQX/20.3.13_40uM_CNQX.mat',4); % done. not very stable (firing rate declines steadily in first 4 hours and continues later)
    sessionConfigs.cnqx40f = newConfig('40uM-CNQX/24.3.13_40uM_CNQX.mat',4); % done
    sessionConfigs.cnqx40g = newConfig('40uM-CNQX/25.9.13_40uM_CNQX.mat',4); % done. this is very different behavior from 24.3.13_40uM_CNQX. Is it really the same condition?
1;
    sessionConfigs.cntx2a = newConfig('2uM-cntx/11.4.13_2uM_cntx.mat',3);
    sessionConfigs.cntx2b = newConfig('2uM-cntx/27.6.13_2uM_cntx.mat',4);    
    sessionConfigs.cntx2c = newConfig('2uM-cntx/30.6.13_2uM_cntx.mat',4);    
    sessionConfigs.cntx2d = newConfig('2uM-cntx/9.7.13_2uM_cntx.mat',4); % first hour are blank due to technical issues
    %sessionConfigs.cntx2e = newConfig('misc/9.7.13_2uM_cntx.mat',4); % PROBLEMATIC. not done. has about 200 negative spike times. Edden - please check this.

    sessionConfigs.aga200a = newConfig('200nM-aga/21.7.13_200nM_aga.mat',4);
    sessionConfigs.aga200b = newConfig('200nM-aga/11.5.14_200nM_aga.mat',3);
    sessionConfigs.aga200c = newConfig('200nM-aga/18.5.14_200nM_aga.mat',3);
     
     
    sessionConfigs.control1 = newConfig('control/27.9.13_control.mat',4); % done
    sessionConfigs.control2 = newConfig('control/1.10.13_15hr control_Mg_HFS.mat',4); % done
    sessionConfigs.control3 = newConfig('control/31.3.13_control.mat',4); % done
    sessionConfigs.control4 = newConfig('control/15.10.13_40hr_cntrl_Abeta_dimer.mat',4); % done
    sessionConfigs.control5 = newConfig('control/17.10.13_control.mat',4); % done

    
    sessionConfigs.b_KO1 = newConfig('1b_KO/7.3.13 1bKO_10uM_bac.mat',4); % done
    sessionConfigs.b_KO2 = newConfig('1b_KO/17.4.13 1bKO_10uM_bac.mat',4); % done
    sessionConfigs.b_KO3 = newConfig('1b_KO/21.4.13 1bKO_10uM_bac.mat',4); % done

    sessionConfigs.a_KO1 = newConfig('1a_KO/20.6.13 1aKO_10uM_bac.mat',4); % done
    sessionConfigs.a_KO2 = newConfig('1a_KO/23.6.13 1aKO_10uM_bac.mat',3); % done
    
    sessionConfigs.more_bac10a = newConfig('more_10uM-bac/16.3.13.mat',4); % done
    sessionConfigs.more_bac10b = newConfig('more_10uM-bac/16.10.12.mat',3); % done
    sessionConfigs.more_bac10c = newConfig('more_10uM-bac/4.12.12.mat',3); % done
    sessionConfigs.more_bac10d = newConfig('more_10uM-bac/9.10.12.mat',3); % done
    sessionConfigs.more_bac10e = newConfig('more_10uM-bac/3.10.13.mat',4); % done
    sessionConfigs.more_bac10g = newConfig('more_10uM-bac/23.11.13.mat',4); % done
    sessionConfigs.more_bac10h = newConfig('more_10uM-bac/20.1.14_10uM_incubated_bac.mat',4); % done

    sessionConfigs.more_control1 = newConfig('more_control/27.9.13_control.mat',4); % done
    sessionConfigs.more_control3 = newConfig('more_control/31.3.13_control.mat',4); % done
    sessionConfigs.more_control4 = newConfig('more_control/15.10.13_40hr_cntrl_Abeta_dimer.mat',4); % done
    sessionConfigs.more_control5 = newConfig('more_control/17.10.13_control.mat',4); % done
    sessionConfigs.more_control6 = newConfig('more_control/28.11.13_control.mat',4); % done
    sessionConfigs.more_control7 = newConfig('more_control/16.1.14_control.mat',4); % done

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% New experiments - Boaz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%       controls with new system 
        sessionConfigs.WTa = newConfig('2015.10.1_WTb_Merged_spikes.mat',3);%this was done with 2X60 runing at the same time

%       MEA Stability (only control hours)from Edden's experiments: 

        sessionConfigs.Cont_SEP_2013 = newConfig('MEA Stability/2013.09.11_Only_cont._Spikes.mat',3);
        sessionConfigs.Cont_OCT_2013 = newConfig('MEA Stability/2013.10.13_Only_cont._Spikes_Sorted.mat',1);
        sessionConfigs.Cont_NOV_2013 = newConfig('MEA Stability/2013.11.12_Only_cont_Spikes.mat',3);

%       48h 10uM CGP incubation tesets:   

        sessionConfigs.CGP10a = newConfig('CGP_Boaz/2014.1.26 CGP/CGPfinal26.1.14(or23).mat',4);
        sessionConfigs.CGP10c = newConfig('CGP_Boaz/2014.2.23/23.2.14CGP10uM_Spikes_Sorted.mat',3);
        sessionConfigs.CGP10d = newConfig('CGP_Boaz/2014.3.14_CGP_10um_Spikes.mat',3);
        sessionConfigs.CGP10d2 = newConfig('CGP_Boaz/14.3.14_10uM_CGP_B.mat',3); %Eddn's analysis
        sessionConfigs.CGP10a2 = newConfig('CGP_Boaz/26.1.14_10uM_CGP_B.mat',3); %Eddn's analysis
        sessionConfigs.CGP10x = newConfig('CGP_Boaz/CGPfinal.mat',3);
%       48h 10uM TeTx incubation tesets:   

        sessionConfigs.TeTx10a = newConfig('2014.4.10_10nM_TeTx_Spikes.mat',4);
        sessionConfigs.TeTx10b = newConfig('2014.5.21 TeTX10nM.mat',4); 
        
%       48h 10uM Baclofen

        sessionConfigs.Bac10a = newConfig('10uMBac_MERGED_Spikes_Sorted.mat',3);
        
%       48h 0.2U chABC

        sessionConfigs.chABCa = newConfig('2014.6.17_chABC0.2U_merged_Spikes.mat',3);
        sessionConfigs.chABCb = newConfig('2015.7.29_b_WTX4_chABC0.2uX24_bac10uM_X48_Spikes_sorted_FFF.mat',4);
        % 0.1-1uM CNQX 
%        
        sessionConfigs.CNQX1a = newConfig('2014.6.26_0.1uMCNQX_Merged_spikes.mat',3);
        sessionConfigs.CNQX1a1 = newConfig('2014.6.26_0.1uMCNQX_Merged_spikes_natural.mat',3);
        sessionConfigs.CNQX1b = newConfig('0.07uM_CNQX_Spikes.mat',3);
        sessionConfigs.CNQX1c = newConfig('0.07uM_CNQX_Spikes-what isthis.mat',3);
        sessionConfigs.CNQX1d = newConfig('1uM_CNQX_1.8.14_Spikes.mat',3);
        
        % 10uM Bicc. 
        sessionConfigs.Bicc10a = newConfig('10uMBicc_1.8_Spikes.mat',3);
        sessionConfigs.Bicc10a1 = newConfig('10uMBicc_1.8_Spikes_natural.mat',3);
        
        % 48h Gabazine 30uM
        sessionConfigs.Gabazine30a = newConfig('Gabazine_30uM_24.8_Merged_spikes.mat',3);
        sessionConfigs.Gabazine30b = newConfig('Gabazine_30uM_2.10_Merged_spikes.mat',3);
        sessionConfigs.Gabazine30bc = newConfig('fdfsfdfsdfsds.mat',3);
        sessionConfigs.Gabazine30c = newConfig('Gabazine_30uM_30.10.14_Merged_UNSORTED.mat',3);
        
        % New sorting for Edden 
        sessionConfigs.XXX = newConfig('.mat',3);
        
        % GBZ + CGP 
         sessionConfigs.GBZCGPa = newConfig('2015.02.17_30uMGBZ_10uMCGP_Spikes.mat',3); % exp of 2015.2.17-19
         sessionConfigs.GBZCGPb = newConfig('2015.03.4_BGZ_CGP_SPIKES_all.mat',3); % exp of 2015.3.2-4 
         sessionConfigs.GBZCGPc = newConfig('2015.03.22 CGP+GBZ Spikes.mat',3); % exp of 2015.03.22-24  
         
        % Teri 
         sessionConfigs.Teri_GBZb = newConfig('2015.04.12-14 Teri+GBZ Spikes.mat',3); 
         sessionConfigs.GBZforTeri = newConfig('2015.04.14-15 30uM GBZ Spikes.mat',3);
         sessionConfigs.GBZforTeric = newConfig('2015.8.10_100uM_Teri_30uM_GBZ.Merged_Spikes0002_sorted.mat',3);
         sessionConfigs.Teri48_Bac = newConfig('2015.1.11_Wt_100uMTeri_10uM_Bac.Mereged.spikes_not_sorted.mat',3);
         sessionConfigs.Teri48_Bac2 = newConfig('2016.4.20_contX4_Teri100uM_ON_Bac10uM48h_Merged_spikes_sorted.mat',4);
         sessionConfigs.TeriON_GBZ48h = newConfig('20016.5.4_100uMTeri_30uMGBZ_Merged_spikes_unsorted.mat',4);
         sessionConfigs.Teri_testTank = newConfig('2016.5.16_X6Cont_Teri100uM(2hj)_spikes_unsorted.mat',3);
         sessionConfigs.Teri_GBZc = newConfig('2016.7.13_2dTeri50uM_incubationX3_GBZ30uM_2hj(3d)_Merged_S0001_unsorted.mat',3);
         sessionConfigs.Teri_A76a = newConfig('2016.7.24_50uMTeri_A7650uM_5uMNN414_300uM1EBIO_TTX_Merged_S_unsorted.mat',4);
         
         %Tau
         sessionConfigs.Tau3ula = newConfig('Tau_spikes_18.6.2015.mat',3);
         
         % DREADD Gq 
         sessionConfigs.Gq_acute = newConfig('15.7.2015_DREADDGq_CNO_Spikes.mat',3);
         sessionConfigs.PVGqa = newConfig('2015.7.29_PV_Gq_4cont_6CNO1uM_45CNO5uM_Spikes_sorted.mat',4);
         sessionConfigs.PVGqb = newConfig('2015.8.12.PV_GqX4_1uM_CNOX12_plus_5uM_CNOX36_Merged_spikes_soted.mat',4);
        
         % DREADD Gi 
         sessionConfigs.PVGia = newConfig('2015.8.24_PV_Gi_cont.1uM_CNO_4cont_49CNO1uM.Merged_spikes_sorted.mat',4);
         sessionConfigs.PVGib = newConfig('2015.10.4_PV_Gi_5uM_CNO_merged_spikes_sorted.mat',4);
         sessionConfigs.PVGic = newConfig('2015.10.4_#2_PV_Gi_5uM_CNO_Merged(3and48)_spikes_sorted.mat',3);
         
         % Channel rhodopsin and SSFO 
         sessionConfigs.PVCha = newConfig('2015.8.30_120_PV_Cre_1ms@20Hz(Cardin2010)_spikes_sorted.mat',3);
         sessionConfigs.PVChb = newConfig('2015.8.30_PV_Chr2_10ms@20Hz_spikes.mat',3);
        
         % M-Current blocker 
         sessionConfigs.M1 = newConfig('2015.10.14_m_curretBLOCkER_contX3_MblockerX49_Merged_spikes_sorted.mat',3);
         sessionConfigs.m1 = newConfig('2015.10.14_m_curret_contX6_BLOCkER_10uMX9_Spikes_sorted.mat',6); %this is full 2h before and full 3h after XE991. it is to test acute effect  
       
         % Lveteriacetame from MAX
         sessionConfigs.Lev1 = newConfig('2015.10.6.WT.LEV35uM.spikes.cutofff.1ms1ms1ms.mat',3);
         
         % Glibenclimde 
         sessionConfigs.Glib1 = newConfig('2015.10.22_WT_10uMGliben_3cont.49Gli_Merged_spikes_sorted.mat',2);
         sessionConfigs.Glib1a = newConfig('2015.10.22_Gliben_full_last_full_2acute_Merged_spikes.mat',3);
         
         % ChR on PV-Cre
         sessionConfigs.ChRa = newConfig('2011.11.7_PV_ChR_Spikes_sorted.mat',3);
         
        % K252 for Ira
        sessionConfigs.K252_a = newConfig('2015.11.10_WT.contX3_200nM_K252a.(24hX1lots2h).Spikes_notSorted.mat',3); 
        sessionConfigs.K252_Bac_b = newConfig('2015.12.10.1 k252a 200nM. Bacl 10uM.merged.Spikes_Merged.mat',3);
        sessionConfigs.K252_Bac_c = newConfig('2016.4.17_200nM_K252a_10uMBac_Merged_spikes_unsorted.mat',4);
        
        % Etomotxir 
        sessionConfigs.Etx_a = newConfig('2015.11.17_contX3_40uM_etomoxirX24X1h_lotsX2h.spikes_notSorted.mat',3);
        sessionConfigs.Etx_b = newConfig('2016.8.29_contX4.40uMEtomoxir(2hj)_Merged_spikes.mat',4); % unsorted
        sessionConfigs.Etx_c = newConfig('2016.9.1contX5_Etx2hj_Bac2jh_Merged_spikes_unsorted.mat',5);
        sessionConfigs.Etx_d = newConfig('2016.9.1ContX4_40uMEtomoxir4hj_10Bac(Etox40uM)4hj.Merged_spikes.mat',4);
        
        % Cont. exp. on new rig
        sessionConfigs.NewCont1 = newConfig('2015.1.7_WT.cont_Merged(2h_jumps)_Spikes_sorted.mat',3);
        
        % KL 001 
        sessionConfigs.KL20 = newConfig('2016.1.12_Wt.8uMKL001_10uMBac.20hMERGE_spikes.mat',3);
        sessionConfigs.KLa = newConfig('2016.1.12_Wt.8uMKL001_10uMBac._washout_Merged_Spikes_unsorted.mat',3);
        sessionConfigs.KLb = newConfig('2016.1.14_WTbX3_8uMKL001X3_GBZ30uMX48h2_Merged_spikes_unsorted.mat',3);  
        
        %cont test - new anti-condensation 2016.2.10
        sessionConfigs.HLa = newConfig('2016.2.9_WT_10umBac_8uM_LH_Merged_spikes_sorted.mat',3);
        sessionConfigs.HLh = newConfig('2016.2.9_WT_10umBac_8uM_LH_spikes_unsorted.mat',2);
        sessionConfigs.HLfull = newConfig('2016.2.9_WT_10umBac_8uM_LH_2hSipks_spikes_unsorted.mat',4);
        sessionConfigs.HL8th9th = newConfig('2016.2.9_WT_10umBac_8uM_LH_8to9hpostHL_only_spikes_unsorted.mat',1);
        sessionConfigs.newConta = newConfig('2016.2.17to22_WTb_cont.only(3hjumps))_spikes_unsorted.mat',2);
        
        % Control Bac 
         sessionConfigs.BacNEWa = newConfig('2016.2.9b_wt_10umBac_Merged(2hskips)_spikes_unsorted.mat',3);
         sessionConfigs.BacNEWb = newConfig('2016.2.25_WTa_10uMBac(3cont.Bac2hJumps)_spikes_unsorted.mat',3);
         sessionConfigs.BacNEWc = newConfig('2016.2.25_WTb10uMBac_3cont_2jumpsBac48h_spikes_unsorted.mat',3);
         sessionConfigs.BacNEWd = newConfig('2016.2.25_WTc10uMBac_3cont._2jumpsBAC48_spikes_unsorted.mat',3);
         
        % PER1/2 Null 
        sessionConfigs.PERa = newConfig('2016.3.9_PER_null_10uMBac_Merged_spikes4.6SDall_unsorted.mat',3); 
        
        % New rig cntrol (+ uridine) 
        sessionConfigs.ContUri = newConfig('2016.4.12_100uM_cont48_Uridine3days_Merged_spikes_sorted.mat',3);
        sessionConfigs.Cont_tankA = newConfig('2016.5.1_cont_WT_Co2tank(noMixer)_spikes_sorted_.mat',3); 
       
        % Antimycin A 
         sessionConfigs.AAa = newConfig('2016.4.23_contX4_1uM_AA_wash_0.5umAA_Merged_spikes_unsorted.mat',4);
         sessionConfigs.AAb = newConfig('2016.8.17a_AA0.5uM.washto0.25uM_10uMBac._X6washout_Spikes_unsorted.mat',4);
        
        % A76 AMPK activator 
        sessionConfigs.A76_a = newConfig('2016.4.27_3Xcont_50uMA76_AMPK_acti_Merged_spikes_sorted.mat',3);
        sessionConfigs.A76_bac_a = newConfig('2016.6.26a_50uMA76_10uMBac_washout_Merged_spikes_unsorted.mat',4);
        sessionConfigs.A76_bac_b = newConfig('2016.8.10b_50uMA76_10uMBac.Merged_spikes_unsorted.mat',3);
        sessionConfigs.A76_bac_c = newConfig('2016.9.8b_contX4_A7650uM(2hj)_10uMBac(2hj)_Merged_S_unsorted.mat',4); 
        sessionConfigs.A76_Teri_a = newConfig('2016.7.13_contX3_50uMA76_Teri50uM(gap)_TTX_Merged_Spikes0001_unsorted.mat',3);
        sessionConfigs.A76_Teri_b = newConfig('2016.8.15b.50uMA76_50uMTeri_Merged_spikes_unsorted.mat',3);
        sessionConfigs.A76_Teri_b = newConfig('2016.8.15b.50uMA76_50uMTeri_Merged_spikes_unsorted.mat',3);
        sessionConfigs.A76_Teri_c = newConfig('2016.9.8_contX4_A7650uM(2hj)_Teri50uM(2hj)_spikes_unsorted.mat',4);
         
        % Uiridine + Teri 
        sessionConfigs.Uri_TeriA = newConfig('2016.5.8_100uMUridine_100uMTeri_Merged_spikes_sorted.mat',3);
        sessionConfigs.Uri_TeriB = newConfig('2016.5.12_100uMUri_100uMTeri_1000uMUri(2nd)_Merged_spikes_unsorted.mat',3);
        sessionConfigs.Uri_TeriBa = newConfig('2016.5.12_100uMUri_100uMTeri_1000uMUri(2nd)_Merged_spikes_sorted.mat',3);
        sessionConfigs.Uri_TeriBb = newConfig('2016.9.23.UMP100uM_Teri50uM_Merged_s_unsorted.mat',3);
        sessionConfigs.Uri_TeriBc = newConfig('2016.9.23b.UMP100uM_Teri50uM_Merged_spikes_unsorted.mat',3);

        % muscimol 
        sessionConfigs.Musci_a = newConfig('2016.5.10_10_4Xcont._uMMus(1hj)_spikes_sorted_active.mat',3);
        sessionConfigs.Musci_b = newConfig('2016.6.26b_10uMMuscimol_washout_Merged(withGap)_spikes_unsorted.mat',3);
        sessionConfigs.Musci_c = newConfig('2016.8.17b_4Xcont_5uMMuscimole.Merged(4hj)_Spikes_unsorted.mat',4);
        
        % shDHODH  
        sessionConfigs.shDHODH_Teri_a = newConfig('2016.5.25_shDHODh_contX3_100uMTeri(1hj)_spikes_unsorted.mat',3);
        sessionConfigs.shDHODH_Teri_b = newConfig('2016.5.28_shDHODH_3Xcont_100uMTeri_24h(1hj)_48h(2hj)_unsorted.mat',3);
        sessionConfigs.shDHODH_Teri_c = newConfig('Data_unsorted.mat',3); %2016.6.2-4
        sessionConfigs.shDHODH_Teri_d = newConfig('2016.7.31_shDHO_teri50um_24h1jh_24h2hj_Merged_spikes_unsorted.mat',3); 
        sessionConfigs.shDHODH_Teri_titr1 = newConfig('2016.6.8_shDHODH_titration_12X5min_Merged_spikes0002_unsorted.mat',3);
        sessionConfigs.shDHODH_Teri_titr2 = newConfig('2016.6.12_shDHODH_titration_Merged_spikes_sorted.mat',3);
        sessionConfigs.shScr_Teri = newConfig('2016.7.28_shScrX3_50uMTeri(1hjX24h2hjX24)_spikes_unsorted.mat',3);
        
        % APP/PS1
        sessionConfigs.APPPS1_GBZa = newConfig('2016.5.30_APPPS1_12Xcont_30uMGBZ24(1hj)48(2hj)_unsorted.mat',3);
        sessionConfigs.APPPS1_GBZb = newConfig('2016.6.5_APPPS1_3Xcont_30uMGBZ_(1hj)_Merged_unsorted.mat',3);
    
        % Teri tirtation on WT 
        sessionConfigs.titTeria = newConfig('2016.6.15b_10_25_50uMTeri(1h each 4h after drug)_Spikes_unsorted.mat',2);
        sessionConfigs.titTeria_full = newConfig('2016.6.15b_Merged_spikes_unsorted.mat',3);
        sessionConfigs.titTeri_b = newConfig('2016.6.21_3Xcont_to50_to100uM_Teri_Merged_spikes_unsorted.mat',3);
        sessionConfigs.titTeri_b2 = newConfig('2016.6.21_cont1h(full)_1h(full)teri_spikes_un.mat',3);
       
        % Trkb IgG
          sessionConfigs.TrkbBac_Aa = newConfig('2016.6.15a_TrkB_IgG5ug10uMBac5ugTrKb_IgGMerged1_spikes_unsorted.mat',3); 
          sessionConfigs.TrkbBac_Ab = newConfig('2016.6.15a_10uMBac_with_5ugTrKb_IgG_again_merged_spikes_unsorted.mat',4); 
          sessionConfigs.TrkbBac_Ac = newConfig('2016.8.10a_igGTrkb10ug_10uMBac.Merged_Spikes_unsorted.mat',3); 

        % STO
         sessionConfigs.STO_GBZa = newConfig('2017.7.18_contX4_GBZ30uM_STO6a09_3uM_2hj_Merged_spikes_unsorted.mat',4);
         
        % TNFa 
         sessionConfigs.TNFa_a = newConfig('2016.8.22b_TNFa_bacl10uM_Merged_spikes_unsorted.mat',4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MAXIM's%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Experiments%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Control
        sessionConfigs.Cont73 = newConfig('12.1.17_Control73h_Spikes.mat',3);
        % AG1024
        sessionConfigs.AG_perliminary = newConfig('2016.10.20_2uMAG.spikes.mat',3);
        sessionConfigs.AG_1uM_Bac10uM_TTX = newConfig('2016.11.27_AG1024_1uM(47)_Bac(48)_TTX_Spikes.mat',3); 
        sessionConfigs.AG_1uM_Bac10uM_Wash = newConfig('2017.1.5_AG+Bac+Wash_spikes.mat',4);
        sessionConfigs.AG_1uM_Bac10uM_Wash2 = newConfig('AG+Bac+Wash_spikes(1ms_1ms_1ms_detection)',18); 
        sessionConfigs.AG_1uM_Bac15uM_Wash = newConfig('2017.3.7_Cont(3)_AG(4)_Bac(48)_Wash(5)_spikes.mat',3); 
        sessionConfigs.AG_1uM_GBZ300uM = newConfig('2017.3.26_Cont(18)_AG(5)_(GBZ)(49)_Spikes.mat',3); 
        sessionConfigs.AG_1uM_GBZ30uM = newConfig('2017.4.5_AG1uM_GBZ30uM_spikes.mat',3); 
        sessionConfigs.AG_1uM_GBZ30uM_unsorted = newConfig('2017.4.3_AG1uM_GBZ30uM_spikes_unsorted.mat',3); 
        sessionConfigs.AG_03uM_in_MEM_unsorted = newConfig('2017.6.21 MEM+AG1024(0.3uM)_Spikes0006.mat',4);
        sessionConfigs.AG_03uM_unsorted = newConfig('2017.8.10_AG0.3uM_Spikes.mat',4);
        sessionConfigs.AG_03uM_Bac_unsorted = newConfig('2017.9.1 AG0.3uM+Bac+Wash+Bac_Spikes.mat',4);
        sessionConfigs.AG_03uM_Bac_unsorted2 = newConfig('2017.9.13_AG0.3uM+Bac10uM+Wash+Bac10uM_Spikes.mat',4);
        sessionConfigs.AG_03uM_Bac_unsorted3 = newConfig('2017.10.8 AG0.3uM+Bac+Wash+Bac_Spikes.mat',4);
        sessionConfigs.AG_03uM_Wash_Bac_unsorted = newConfig('2017.10.30_AG0.3uM+Wash+Bac_Spikes.mat',4);
         
        %GBZ
        
        sessionConfigs.Cont_GBZ30uM = newConfig('2017.3.7_Cont(21)_GBZ30uM(48)_Spikes.mat',3); 
        %Bac
      sessionConfigs.BacCont1 = newConfig('Bac10uM(48).mat',3);
      sessionConfigs.BacCont2 = newConfig('Bac10uM_spikes001.mat',16);
      sessionConfigs.BacCont3 = newConfig('2017.2.1_Bac10uM(48h)_spikes.mat',3);
      sessionConfigs.BacCont4 = newConfig('2017.2.8 Bac20uM_Bath_TTX.mat',4);
       sessionConfigs.BacCont5 = newConfig('2017.2.23_Cont(3)_Bac10uM(47)_spikes.mat',3);
      sessionConfigs.BacIncubator = newConfig('2017.2.22_Bac(incubator)_Spikes.mat',3);
        %shIGF1R
         sessionConfigs.KD_Bac1 = newConfig('2017.5.4_shIGF1R+Bac+Wash_Spikes0001.mat',3);
         sessionConfigs.KD_Bac2 = newConfig('2017.6.4_shIGF1R+Bac(66)+Wash(16)_Spikes0002.mat',3);
         sessionConfigs.KD_Bac3 = newConfig('2017.6.29_shIGF1R+Bac_Spikes.mat',3);
         sessionConfigs.KD_Bac4 = newConfig('2017.7.5 IGF1R-KD+Bac+Wash_Spikes.mat',3);
        %shSCR
         sessionConfigs.SCRBac1 = newConfig('2017.6.29_Scr+Bac_Spikes.mat',3);
         %Valproate
         sessionConfigs.VPA1 = newConfig('2017.5.15_Sodium_valproate_Spikes.mat',3);
         sessionConfigs.VPA2 = newConfig('2018.1.14_Valproate_600uM_Washout_Spikes.mat',4);
         %Retigabine
         sessionConfigs.Reti1 = newConfig('2017.5.22_Retigabine10uM.mat',3);
         sessionConfigs.Reti2 = newConfig('2017.7.16_Retigabine10uM_Spikes.mat',4);
         
         %ENBA
         sessionConfigs.ENBA1 = newConfig('2017.8.6_(4hCont)ENBA10nM (2X2)50nM(2X24)_Spikes.mat',4);
         sessionConfigs.ENBA2 = newConfig('2017.9.1_ENBA10nM_Spikes.mat',4);
         sessionConfigs.ENBA_GBZ = newConfig('2017.12.5_ENBA10nM_GBZ10uM_Spikes.mat',3);
         sessionConfigs.ENBA1_5_GBZ10 = newConfig('2017.12.20ENBA_1nM_5nM+GBZ10uM_spikes.mat',4);
         %Muscarine
           sessionConfigs.Musc1 = newConfig('2017.7.24_Muscarine10uM_Spikes.mat',4);
           
         %XE991
          sessionConfigs.XE = newConfig('2017.11.9_XE991_10uM_Spikes.mat',4);
          sessionConfigs.XEboaz = newConfig('2017.7.20A_XE199_Merged_spikes_un.mat',3);
           sessionConfigs.XEsorted1 = newConfig('2017.7.20_XE991_Spikes.mat',3);
            sessionConfigs.XEsorted2 = newConfig('2017.7.24_XE991_Spikes.mat',3);
             sessionConfigs.XEsorted3 = newConfig('2017.8.25_XE991_Spikes.mat',3);
              sessionConfigs.XEsorted4 = newConfig('2017.11.9_XE991_10uM_1h_Spikes.mat',3);
           
         %Caffeine
            sessionConfigs.Cafe = newConfig('2017.11.6_Caffeine 10uM _Spikes.mat',4);
         %A76
          sessionConfigs.A76GBZ = newConfig('2017.10.30_A7650uM_GBZ10uM_Spikes.mat',4);
          
         %Teriflunomide
         sessionConfigs.TeriTest = newConfig('2018.1.7_Teri50_100Test.mat',1);
         
         %4AP tests (Boaz)
          sessionConfigs.the4APa = newConfig('2018.1.18bContX4(1hj))_50uM4AP(1hj_GAP_X1hj)_MergedS_un.mat',4);
          
          %Neta in-vivo
            sessionConfigs.Neta1 = newConfig('1_5_iso_begfore_ptz_4stdev-1.mat',1);
         sessionConfigs.Neta2 = newConfig('1_5_iso_begfore_ptz_4stdev-2.mat',1);
         sessionConfigs.Neta3 = newConfig('Block-5.mat',1);
        
         %Anto MEAs
         
          %APPPS1 FOR DANIEL PAPER
         sessionConfigs.APPPS1_GBZ1=newConfig('2019.2.27_Anto_A_APPPS_GBZ_spike sorted.mat',3);
         sessionConfigs.APPPS1_GBZ2=newConfig('2019.1.25_APPPS1_48hGBZ_spike sorted.mat',4);
         sessionConfigs.APPPS1_GBZ3=newConfig('Spike sorted_Exp57c_2h_new.mat',4);
         sessionConfigs.APPPS1_GBZ4=newConfig('25072018_APPPS1_5hbsln-48hGBZ.mat',5);
         %sessionConfigs.APPPS1_GBZ5=newConfig('Exp57b_APPPS1_4hBSL_every2h.mat',4);
         sessionConfigs.APPPS1_GBZ6=newConfig('2019.3.27_APPPS1(6h)_30uMGBZ_spike sorted0002.mat',5);
         sessionConfigs.Cntrl_GBZ1=newConfig('Gabazine_30uM_2.10_Merged_spikes.mat',3);
         sessionConfigs.Cntrl_GBZ2=newConfig('Gabazine_30uM_24.8_Merged_spikes.mat',3);
         sessionConfigs.Cntrl_GBZ3=newConfig('Gabazine_30uM_30.10.14_Merged_UNSORTED.mat',3);
         %sessionConfigs.Cntrl_GBZ4=newConfig('2015.02.17_30uMGBZ_10uMCGP_Spikes.mat',3);
         %sessionConfigs.Cntrl_GBZ5=newConfig('2017.7.18_contX4_GBZ30uM_STO6a09_3uM_2hj_Merged_spikes_unsorted.mat',3);
         sessionConfigs.Cntrl_GBZ6=newConfig('2015.04.14-15 30uM GBZ Spikes.mat',3);
         %sessionConfigs.DIAZ_GBZ=newConfig('2020.4.29_A_5uM DIAZ_5uM GBZ_spikes.mat',4);
         %sessionConfigs.Cntrl_GBZ_HD8=newConfig('2014.12.08_4hBSL_gbz_spike sorted.mat',4);
         %sessionConfigs.Cntrl_GBZ_HD8_b=newConfig('2014.12.08_4hBSL_gbz_spikes.mat',4);
         sessionConfigs.Cntrl_GBZ_Max=newConfig('2021.9.17_NoCre_GBZ10uM_Spikes.mat',4);
         sessionConfigs.APP_GBZ1_SS=newConfig('2019.2.27_APPPS_GBZ1_NEW_spike sorted.mat',3);
         sessionConfigs.APP_GBZ2_SS=newConfig('2019.1.25_APPPS1_GBZ2_NEW_spike sorted.mat',4);
         %sessionConfigs.APP_GBZ3_SS=newConfig('57C__APPPS_GBZ3_spike sorted',3);
         sessionConfigs.APP_GBZ6_SS=newConfig('2019.3.27_APPPS1(6h)_GBZ6_NEW_spike sorted.mat',5);
         sessionConfigs.Cntrl_A2_SS=newConfig('2021-09-23T17-05-04_mea4_12hBSL_gbz30uM_spikes SORTED.mat',12);
         sessionConfigs.Cntrl_A3_SS=newConfig('2021-09-23T17-08-38_mea1_12hBsl_30uM GBZ_spikes SORTED.mat',12);
         sessionConfigs.Cntrl_GBZ_Max_SS2=newConfig('Max_GBZ2.mat',4);
         sessionConfigs.Cntrl_A2_SS2=newConfig('Cntrl_mea4_spikes.mat',12);
         sessionConfigs.Cntrl_A3_SS2=newConfig('Cntrl_mea1_spikes.mat',12);
         
         %
         sessionConfigs.Exp50_Analysis = newConfig('Exp50_3hBsl_1hCompleteARC_2hARC_1hUK5099_1hBac_sorted.mat',3);
        sessionConfigs.Exp53_Analysis = newConfig('2018.5.1_spike soreted_analised.mat',4);
        sessionConfigs.Exp53b_Analysis = newConfig('2018.5.1_spike sorted_5.5SD.mat',4);
        sessionConfigs.Exp55_Analysis= newConfig('Exp55_spike sorted_4h_bsln+uk5099_5h1_17h10_7h50_17h100.mat',4);
         sessionConfigs.Exp56a_Analysis= newConfig('Exp56_APPPS1_4hbsln_30uMGBZ_every2h.mat',4);
         sessionConfigs.Exp56b_Analysis= newConfig('Exp57b_APPPS1_4hBSL_every2h.mat',4);
         sessionConfigs.Exp56c_Analysis= newConfig('Spike sorted_Exp57c_2h_new.mat',4);
         sessionConfigs.Exp57b2_Analysis= newConfig('Exp57b_sorted_2h(4)',4);
         sessionConfigs.Exp58_DAB=newConfig('Exp58_spike sorted_Channel a',4);
         sessionConfigs.Exp58b_DAB=newConfig('Exp58b_spike sorted',4);
        sessionConfigs.Exp58c_Analysis=newConfig('Exp58c',4);
        sessionConfigs.Exp60_b_STO= newConfig('2018.7.5_Analysis_4hsto_48hgbz_30uM.mat',4);
        sessionConfigs.Exp50_new_Analysis=newConfig('Exp50_new sorting.mat',3);
        sessionConfigs.Exp54_new_Analysis=newConfig('Exp54_spike sorted_4h_bsln+uk5099_5h1_17h10_7h50_17h100.mat',4);
        sessionConfigs.Exp61_APPPS1=newConfig('25072018_APPPS1_5hbsln-48hGBZ.mat',5);
        sessionConfigs.shPDZD8_baclofen=newConfig('2018.8.15_shPDZD8_10uMBac.mat',4);
         sessionConfigs.shPDZD8_GBZ=newConfig('2018.8.15_shPDZD8_30uMGBZ.mat',4);
         sessionConfigs.scrmPDZD8=newConfig('2018.8.17_scramble_4hBsln_48hBac0001.mat',4);
         sessionConfigs.Eto_Bac_MEAa=newConfig('2018.Analysis_MeaA_100uM Etomoxir_Baclofen.mat',4);
         sessionConfigs.Eto_Bac_MEAb=newConfig('2018_Analysis_MEAB_etomoxir+Baclofen',4);
         sessionConfigs.shCPT1c_2ul=newConfig('2018.09.04_Analysis_shCPT1c_2ul_4hbsln_baclofen.mat',4);
         sessionConfigs.shCPT1c_4ul=newConfig('2018.04.09_Analysis_shCPT1c_4ul_4hbsln_baclofen.mat',4);
         sessionConfigs.shPDZD8_Bac_II=newConfig('2018.09.15_shPDZD8_10uMBac(b).mat',4);
         sessionConfigs.Eto10uM__Bac=newConfig('2018.10.26_10uMEto_10uMBac.mat',4);
         sessionConfigs.Dandrolene_10uM__Bac=newConfig('2018.10.26_spike sorted_5hbln_10uMDantrolene_10uMBac.mat',5);
         sessionConfigs.Eto100uM_GBZ30uM=newConfig('2018.11.6_spike sorted_4hbln_100uMEtomoxir_30uMgbz.mat',4);
          sessionConfigs.Dand_BacII=newConfig('2018.11.19_spike sort_4hbsln_10uMDand(1h)25h_10uMBac(2h).mat',4);
          sessionConfigs.Perhexilin100uM=newConfig('2018.7.5_spike sorted_4hsto_100uM Perhexilin.mat',4);
          sessionConfigs.shNCLX_GBZ=newConfig('2018.11.29_10uL shNclx_30uM GBZ.mat',4);
          sessionConfigs.shNCLX_GBZ_II=newConfig('2018.12.9_14uL shNclx_30uM GBZ.MAT',4);
          sessionConfigs.Control_GBZ=newConfig('9.12.18_Control_30uMGBZ(2h)_spike sorted.mat',4);
          sessionConfigs.Ketamine=newConfig('2019.1.1_spike sort_4hbsln_50uMKetamine(2h)24h_10uMBac(2h).mat',4);
         sessionConfigs.Ketamine_2=newConfig('2019.1.7_Anto_50uM ketamine_10uM baclofen_spike sorted.mat',4);
         sessionConfigs.Eto40uM__Bac=newConfig('2019.1.1_spike sort_4hbsln_50uMEtomoxir(2h)24h_10uMBac(2h).mat',4);
         sessionConfigs.Ketamine_3=newConfig('2019.1.22_spike sort_4hbsln_50uMKetamine(2h)24h_10uMBac(2h).mat',4);
         sessionConfigs.APPPS1_baclofen=newConfig('2019.1.25_APPPS1_48hBaclofen_spike sorted.mat',4);
         sessionConfigs.APPPS1_GBZ4=newConfig('2019.1.25_APPPS1_48hGBZ_spike sorted.mat',4);
         sessionConfigs.APPPS1_bac2=newConfig('2019.2.4_APPPS1_72hBaclofen_spike sorted.mat',4);
         sessionConfigs.Eto40uM_Bac_050219=newConfig('2019.2.5_Anto_40uM etomoxir_10uM bac.mat',4);
         sessionConfigs.KETAMINE_4_Bac=newConfig('2019.2.5_Anto_50uM KetA(6H)_10uM bac_SPIKE SORT.mat',4);
         sessionConfigs.Atest=newConfig('Anto_testA',2);
         sessionConfigs.shCPT1c_3=newConfig('2019.2.10_shCPT1c_10uMBac_spike sorted.mat',10);
         sessionConfigs.DS_TBOA1=newConfig('18.02.19_Dravet_6hbsln_TBOA_spike sorted.mat',6);
         sessionConfigs.DS_TBOA2=newConfig('17.02.19_B_Dravet_6hbsln_TBOA_spike sorted.mat',6);
         sessionConfigs.APPPS1_GBZ5=newConfig('2019.2.27_Anto_A_APPPS_GBZ_spike sorted.mat',4);
         sessionConfigs.APPPS1_ketamine_BAC=newConfig('2019.3.2_Anto_A_APPPS150uM KETAMINE_BACLOFEN_SPIKE SORTED.mat',4);
         sessionConfigs.APPPS1_ketamine6h=newConfig('2019.2.5_Anto_50uM KetA(6H)_10uM bac_SPIKE SORT.mat',4);
         sessionConfigs.Perhexiline_1uM=newConfig('2019.2.6_Anto_1uM perhexiline_spike sorted.mat',4);
         sessionConfigs.Perhexiline_1uM_bac=newConfig('2019.2.6_Anto_1uM perhexiline_bac_spike sorted.mat',4);
         sessionConfigs.Propofol_5uM_1h=newConfig('3hbsln_1hprop.mat',3);
         sessionConfigs.Propofol_5uM_Exp1=newConfig('6hbsln_1h_I_5uM_prop_spike sorted.mat',6);
         sessionConfigs.shCPT1c_4_GBZ=newConfig('2019.3.10_Anto_shCPT1c_6hbsln_30uM gbz_spike sort0001.mat',6);
         sessionConfigs.Propofol_10uM=newConfig('Propofol_10uM_spike sorted.mat',6);
         sessionConfigs.Test1=newConfig('test A_spike sort',3);
         sessionConfigs.APP_Ps1_Gbz_cont=newConfig('2019.3.19_APPPS1(5h)_gbz_contamination0001.mat',5);
         sessionConfigs.APP_PS1_baclofen3=newConfig('2019.3.19_APPPS1(6h)_baclofen.mat',6);
         sessionConfigs.APP_PS1_GBZ6=newConfig('2019.3.27_APPPS1(6h)_30uMGBZ_spike sorted0002.mat',6);
         sessionConfigs.shCPT1c_5_GBZ=newConfig('3.4.19_shCPT1c_6hbsln_gbz.mat',6);
         sessionConfigs.Perhexiline_1uM_bac2=newConfig('4.3.19_Perhexiline_baclofen_spike sort.mat',4);
         sessionConfigs.Mg_test=newConfig('3.4.19_mg test.mat',3);
         sessionConfigs.Prop_test=newConfig('3.4.19_prop test0001.mat',3);
         sessionConfigs.Prop_bac=newConfig('2019.4.4_4hbsln_10uM_prop_10uM bac_spike sorted.mat',4);
         sessionConfigs.Mg_test2=newConfig('2019.4.14_Anto_b_1.2mM_Mgspike sort.mat',3);
         sessionConfigs.Mg_bac=newConfig('2019.4.4_6hbsln_1uM_MgCl2_10uM bac_spike sorted0001.mat',6);
         sessionConfigs.Diazepam10uM=newConfig('2019.4.14_Anto_A_10uM diazepam_bac_merged_spike_sort.mat',6);
         sessionConfigs.Mg_bac2=newConfig('2019.4.14_6hbsln_1uM_MgCl2_10uM bac_spike sorted.mat',6);
         sessionConfigs.Control_Bac=newConfig('2019_4_17_Anto_A_control_10Um BAC_spike sort.mat',6);
         sessionConfigs.APPPS1_Bac4=newConfig('2019_4_17_Anto_B_APPPS1_6h_bsln_10uM bac_spike sort.mat',6);
         sessionConfigs.Stability_04_19=newConfig('Stability_1.2mM MgCl2_spike sort.mat',6);
         sessionConfigs.Stability_BAC_04_19=newConfig('2019.4.29_anto_A_6h bsln_10uM bac_spike sort.mat',6);
         sessionConfigs.Diazepam_test=newConfig('2019.5.28(A)_Cont_5uM diazepam_SPIKE SORT_TEST.mat',3);
         sessionConfigs.A76_test=newConfig('2019.5.28(B)_Cont_50uM A76_SPIKE SORT 3+3(2)_TEST.mat',3);
         sessionConfigs.A76_test2=newConfig('2019.5.28(B)_Cont_50uM A76_SS_TEST.mat',1);
         sessionConfigs.Diazepam5uM_Bac=newConfig('2019.5.28(A)_Cont4h_5uM diazepam5h(1h)_10um bac_spike sort.mat',4);
         sessionConfigs.A76_baclofen=newConfig('2019.5.28(B)_Cont4h_50uM A76_6h(1h)_bac_spike sort.mat',4);
         sessionConfigs.isof_test=newConfig('iso_test_ss.mat',2);
         sessionConfigs.isof150_test=newConfig('test150_spikesort.mat',2);
         sessionConfigs.Isoflurane_test=newConfig('5.6.19_isoflurane_test_3hbsl_1%40mLxmin(3h_1h)_1%_150mlmin(3h_1h)_washout(5x1h_2h).mat',3);
         sessionConfigs.Isoflurane_test2=newConfig('2019.6.4_isoflurane_test_spike sorted.mat',4);
         sessionConfigs.KETAMINE4=newConfig('2019.2.5_Anto_50uM KetA(6H)_10uM bac_SPIKE SORTED(REAL).mat',4);
         sessionConfigs.KETAMINE48h=newConfig('2019.6.12_4H BSLN_50uM ket(2H_27FILES)_50uM A76.mat',4);
         sessionConfigs.KETAMINE48h_2=newConfig('2019.6.12_4H BSLN_50uM ket(2H_27FILES)_50uM A76_SS2.mat',4);
         sessionConfigs.Ketamine_TBOA=newConfig('2019.6.12_4HBSLN_8H50uM ketamine_10uM TBOA.mat',4);
         sessionConfigs.APV_Ketamine=newConfig('2019.6.16_Anto_4H BSLN_10uM AP-V(25files_2H)_50uM KET_SS.mat',4);
         sessionConfigs.KETAMINE_GBZ=newConfig('2019.7.10_Anto_4hbsln__5h 50uM ketamine_10uM gbz_ss.mat',4);
         sessionConfigs.KETAMINE_TEST=newConfig('Test_b_ket.mat',2);
         sessionConfigs.Ketamine_testB=newConfig('2019.7.29_B_new_frozen ketamine_4h_bsln_50+50.mat',4);
         sessionConfigs.Ketamine_testA=newConfig('2019.7.29_Anto__A_RT 50uM KET_NEW+50_SPIKE S.mat',4);
         sessionConfigs.Ketamine_APV_B=newConfig('2019.7.29_4h_100uM keta+4h APV_meaB.mat',4);
         sessionConfigs.Ketamine_APV_A=newConfig('2019.7.29_4h_100uM keta+4h APV_mea_A.mat',2);
         sessionConfigs.Diazepam5uM_2=newConfig('2019.8.2_5uM diazepam.mat',4);
         sessionConfigs.Diazepam5uM_bac_2=newConfig('2019.8.2_4hbsln_5uM diazepam_Bac10uM_ss.mat',4);
         sessionConfigs.Ketamine_inMgCl=newConfig('2019.8.17_Anto_A_4H_BLSN_1.2mM MgCl2_50uM KETA_SS.mat',4);
         sessionConfigs.Ketamine_DIV10=newConfig('2019.8.17_Anto_B_DIV10_50uM KETA_SS.mat',4);
         sessionConfigs.Ketamine_DIV10=newConfig('2019.8.17_Anto_B_DIV10_50uM KETA_SS2.mat',4);
         sessionConfigs.Ketamine_TBOA=newConfig('2019.8.17_4h_1.2mM MgCl_28h_50uM KET_48h_5uM TBOA_ss.mat',4);
         sessionConfigs.Ketamine_TBOA2=newConfig('2019.8.17_DIV10_4hBSL_48h_50uM KETA_24h_5uMTBOA_ss.mat',4);
         sessionConfigs.Diazepam5uM_Bac3=newConfig('2019.8.21_Anto_4hblsn_8x5uM Diazepam(2h)_10uMbaclofen_SS.mat',4);
         sessionConfigs.Diazepam5uM_TBOA=newConfig('2019.8.21_Anto_4hBsln_5x50uM diazepam_5uM TBOA_ss.mat',4);
         sessionConfigs.Diazepam5uM_Bac3b=newConfig('2019.8.21_Anto_4hblsn_8x5uM Diazepam(2h)_10uMbaclofen_SS2.mat',4);
         sessionConfigs.Diazepam5uM_TBOA_2=newConfig('26.8.19_5uM Diazepam(pre-inc)6h_5uM TBOA_mat.mat',6);
         sessionConfigs.Diazepam5uM_TBOA_3=newConfig('1.9.19_5uM Diazepam16h_5uM TBOA_ss.mat',4);
         sessionConfigs.Diazepam5uM_TBOA_4=newConfig('1.9.19_B_5uM Diazepam_5uM TBOA_ss.mat',4);
         sessionConfigs.Diazepam5uM_TBOA_5=newConfig('2019.10.10_diaz_TBOA.mat',4);
         sessionConfigs.Ketamine_TBOA3=newConfig('2019.10.6_A_1.2mM MgCl(4h)_50uM keta_5uM TBOA_SS.mat',4);
         sessionConfigs.Diazepam5uM_washout=newConfig('2019.10.15_4h bsln_5uM diazepam(30h circa)_4hwashout(1h)_new.mat',4);
         sessionConfigs.eEF2KO_Bac=newConfig('2019.10.16_eEF2KO_8hbsln_10uMBaclofen.mat',8);
         sessionConfigs.eEF2KO_TBOA=newConfig('2019.10.20__eEF2KO_8HBSLB_10uM TBOA.mat',8);
         sessionConfigs.Isoflurane_TBOA=newConfig('201910.20_4hBSLN__Iso_10uM TBOA(2h)_wash(1h).mat',4);
         sessionConfigs.TNFa_Diaz=newConfig('20101030_TNFa KO_8HBSLN_6Hdiazepam_bac.mat',8);
         sessionConfigs.eEF2KO_Bac2=newConfig('20191030_eEF2KO_8h BSLN_bac.mat',4);
         sessionConfigs.Diazepam5uM_TBOA_6=newConfig('20191103_4HBSLN_5uM DIAZ(4X2h)_GAP TILL 11.15_10uM TBOA.mat',4);
         sessionConfigs.eEF2KO_TBOA2=newConfig('20191103_eEF2KO7HBSLN_10uM TBOA(2H)_WASHOUT(2H).mat',7);
         sessionConfigs.ENBA_DIAZ=newConfig('2019.11.10_4hbsln_10nM ENBA(3h)_diaz(2h)_ss.mat',4);
         sessionConfigs.TERI_DIAZ=newConfig('2019.11.10_4hbsln_50uM teri(4h)_5uM diaz(2h)_ss.mat',4);
         sessionConfigs.Isoflurane48h=newConfig('2019.11.20_5h bsln_48h_0.2xcent_iso35mlxmin_ss0001.mat',5);
         sessionConfigs.TERI_DIAZ2=newConfig('2019.11.25_Anto_A_50uM Teri_5uM DIAZ.mat',4);
         sessionConfigs.ENBA_DIAZ2=newConfig('2019.11.25_10nM_ENBA(6h)_DIAZ.mat',4);
         sessionConfigs.ENBA_DIAZ2b=newConfig('2019.11.25_10nM_ENBA(6h)_DIAZ_diff analysis.mat',4);
         sessionConfigs.TNFa_TBOA=newConfig('2019.12.1_TNFa6h BSLN_(24h)5uM DIAZ_10uM TBOA.mat',6);
         sessionConfigs.TNFa_BAC2=newConfig('2019.12.1_TNFa_6h BSLN_5uM DIAZ(6H)_10uM BAC.mat',6);
         sessionConfigs.Midazolam=newConfig('2019.11.8_6hbsln_3x1uM midazolam.mat',6);
         sessionConfigs.APP_PS1_KETA_BAC=newConfig('2019.12.11_APPPS1_6hbsln_3x1h_ketamine_Bac(2h).mat',6);
        sessionConfigs.eEF2KO_KETA_a=newConfig('2019.12.21_A_eEF2KKO_6h bsln_13h ketam(1h)_Bac_spikes.mat',6);
        sessionConfigs.eEF2KO_KETA_b=newConfig('2019.12.21_B_eEF2KKO_6h bsln_13h ketam(1h)_Bac_spikes.mat',6);
        sessionConfigs.Domitor_1=newConfig('2019.12.18_4hbsln_domitor(2h)_baclofen_spikes.mat',4);
        sessionConfigs.Domitor_5uM_Bac=newConfig('2019.12.30_Anto_B_5uM domitor_bac_spikes.mat',4);
        sessionConfigs.eEF2KO_KETA_c=newConfig('2019.12.30_6hbsl_eEF2KKO_3x50uM ket(2h)_10uM bac_sp.mat',6);
        sessionConfigs.A76_ketamine=newConfig('2020.1.12_50uM A76 in 1.2mM Mg(8h_1h)_50uM keta_s.mat',8);
        sessionConfigs.A76_diaz=newConfig('2020.1.12_50uM A76_diaz.mat',8);
       sessionConfigs.Domitor_5to10uM_TBOA=newConfig('2020.23.1_domitor 5uM_10uM_TBOA.mat',4);
       sessionConfigs.WT_ket_bac=newConfig('2020128_WT_1.2mM MgCl2(6h_1h)_KET50uM(x5_2h)_bac_S.mat',6);
       sessionConfigs.eEF2KKO_4_ket_bac=newConfig('2020128_eEF2K KO_50uM KET(2h-X6h)_bac.mat',6);
       sessionConfigs.eEF2KKO_5_ket=newConfig('2020131_A_1EF2K KO_1.2mM MgCl2 BSLN_6h_50uM KET(2h)_S.mat',6);
       sessionConfigs.eEF2KKO_5_ket_bac=newConfig('2020131_A_1EF2K KO_1.2mM MgCl2 BSLN_6h_50uM KET(2h)_BAC_SPIKES.mat',6);
       sessionConfigs.eEF2KKO_6_ket_bac=newConfig('2020131_B_1EF2K KO_1.2mM MgCl2 BSLN_6h_50uM KET(2h)_BAC_SPIKES.mat',6);
       sessionConfigs.eEF2KKO_1_diaz_bac=newConfig('2020.2.4_e4hbsln_EF2KO_o.n.5uM DIAZ_BAC_SPIKES.mat',4);
       sessionConfigs.ENBA_DIAZ_2=newConfig('2020.2.18_B_4HBSLN_3X100nM ENBA(2h)_5uM DIAZ(2H)_SPIKES.mat',4);
       sessionConfigs.ENBA_keta_1=newConfig('2020.2.20_4hbsln_3xENBA100uM(gap18h)_50uM keta24h_spikes.mat',4);
       sessionConfigs.eEF2KKO_K252=newConfig('2020.2.23_6hbsln_eEF2K KO_200nM K252(24h_2h)_BAC_spikes.mat',6);
       sessionConfigs.Domitor_bac=newConfig('2020.2.18_4hbsln__5uM DOMITOR_BAC3d_spikes.mat',4);
       sessionConfigs.eEF2KKO_BAC3=newConfig('2020.2.23_eEF2K KO_6h bsln_BAC_spikes.mat',6);
       sessionConfigs.TNFa_BAC3=newConfig('2020.2.9_6hbsln_TNFa_5x2h_DIAZ_bac_spikes.mat',6);
       sessionConfigs.KN93_DIAZ=newConfig('2020.3.5_Anto_B_Ctrl_5uM KN-93_Diazepam_spikes.mat',4);
       sessionConfigs.K252_Bac=newConfig('2020.3.1_Anto_4HBSL_24h_K252(2H)_bac_spikes.mat',4);
       sessionConfigs.KN93_Keta=newConfig('2020.3.10_6h bsln_1.2mMMg_(22h_2h)5uM kn93_24h_KETA_SPIKE.mat',6);
       sessionConfigs.STO_Diaz=newConfig('2020.3.10_4HBSLN_3uM STO(2h_22h)__5uM_DIAZ24H_SPIKES.mat',4);
       sessionConfigs.KN92_Bac=newConfig('2020.3.18_4h BSLN_5uM KN92_24h(2h)_2dBAC_SPIKES.mat',4);
       sessionConfigs.KN93_Bac=newConfig('2020.3.18_4hBSLN_12hKN 93 5uM(2h)_2dBAC_SPIKES.mat',4);
       sessionConfigs.Cpt1cKO_Bac_1=newConfig('2020.3.24_CPT1cKO_bac.mat',4);
       sessionConfigs.Cpt1cKO_Bac_2=newConfig('2020.3.27_CPT1cKO_4h BSLN_10uMBAC_2.mat',4);
       sessionConfigs.Cpt1cKO_GBZ_1=newConfig('2020.3.26_CPT1cKO_4h BSLN_30uMGBZ_1.mat',4);
       sessionConfigs.Cpt1cKO_GBZ_2=newConfig('2020.3.30_CPT1C_4h bsln_Gbz10uM_2.mat',4);
        sessionConfigs.Cpt1cKO_Bac_3=newConfig('2020.3.30_CPT1C_A_Bac10uM_SPIKES.mat',4);
        sessionConfigs.TatCN21_Diaz_1=newConfig('2020.4.12_4hbsln_5uMCN21(6x1h_2h)_20uM(7x1h)_Diaz_SPIKES.mat',4);
        sessionConfigs.TatCntrl_Diaz_1=newConfig('2020.4.12_4hbsln_5uMCNTRL(6x1h_2h)_20uM(7x1h)_Diaz_SPIKES.mat',4);
        sessionConfigs.Cpt1cKO_GBZ_3=newConfig('2020.3.26_CPT1C_4h bsln_Gbz30uM_3.mat',4);
        sessionConfigs.KN93_DIAZ2=newConfig('2020.4.15_A_5uM KN93(24h_2h)_5uM DIAZ(24h_2h)_spikes.mat',4);
        sessionConfigs.KN92_DIAZ=newConfig('2020.4.15_B_5uM KN92(2H_24H)_5uM DIAZ(2H_24H).mat',4);
        sessionConfigs.TatCntrl_bac=newConfig('2020-04-23T11-01-2414.4.20_tatCNTRL_BAC.mat',4);
        sessionConfigs.KN93_Keta2=newConfig('2020-04-26T12-24-1523.4.20_KN93_KET.mat',4);
        sessionConfigs.Diaz_GBZ=newConfig('2020.4.29_A_5uM DIAZ_5uM GBZ_spikes.mat',4);
        sessionConfigs.Diaz_GBZ2=newConfig('2020.4.29_B_5uM DIAZ_5uM GBZ_spikes.mat',4);
        sessionConfigs.WT_BAC_newMEA=newConfig('5.5.20_24hbsln_bac.mat',4);
        sessionConfigs.Propofol_Bac_new=newConfig('2020.5.10_4hbsln_5uM PROPOFOL_Bac10uM_SPIKES.mat',4);
        sessionConfigs.Domitor_5uM_Bac3=newConfig('2020.5.10_4hbsln_24h_5uM DOMITOR_Bac10uM_SPIKES.mat',4);
        sessionConfigs.Domitor_5uM_Bac4=newConfig('2020-05-15T11-02-139.5.20_4hbsln_24h(2h)Domitor_Bac.mat',4);
        sessionConfigs.DHODHR135C=newConfig('2020.5.17_DHODH_R135C_4hBSLN__Bac_SPIKES.mat',4);
        sessionConfigs.DHODH_WT=newConfig('2020.5.17_DHODH_WT_4hBSLN__Bac_SPIKES.mat',4);
        sessionConfigs.Diaz_GBZ3=newConfig('2020.5.26_5uM DIAZ(gap 24h_last 2 points_1h)_5uM GBZ_spikes.mat',4);
        sessionConfigs.Propofol_48h=newConfig('2020.5.27_2h_2h(24h gap)baseline_5uMPROPOFOL_48h_spikes.mat',4);
        sessionConfigs.diaz_ket=newConfig('2020.6.2_Anto_B_5uM diaz_50uM keta_spikes.mat',4);
        sessionConfigs.STO_diaz2=newConfig('2020.6.4_Anto_A3uM_STO(gap)_5uM diaz_spikes.mat',4);
        sessionConfigs.Diaz_GBZ4=newConfig('2020-06-10T19-47-44Diaz_GBZ.mat',4);
        sessionConfigs.Keta_Diaz=newConfig('2020.6.11_1.2mMMgCl2_keta_bac_spikes.mat',4);
        sessionConfigs.Keta_Diaz_b=newConfig('2020.6.11_1.2mMMgCl2_keta_bac_spikes2.mat',4);
        sessionConfigs.Prop_3=newConfig('2020-06-12T12-00-482020.6.12_5uM propofol.mat',4);
        sessionConfigs.Keta_Diaz_full=newConfig('2020.6.11_1.2mMMgCl2_keta_bac_spikes full exp.mat',4);
        sessionConfigs.Diaz_Bac4=newConfig('2020.6.11_4h bsl_24h_5uM diaz_bac_spikes.mat',4);
        sessionConfigs.CaMKIIabKD_BAC_1=newConfig('2020.6.16_CamKIIab KD_MEA_A_bac_spikes.mat',6);
        sessionConfigs.CaMKIIabKD_BAC_2=newConfig('2020.6.16_CamKIIab KD_MEA_B_bac_spikes.mat',6);
        sessionConfigs.CaMKIIabKD_DIAZ=newConfig('2020.6.23_Anto_6hbsln_4hgap(1-2)__CaMKIIKD_24hdiaz_spikes.mat',6);
        sessionConfigs.CaMKIIabKD_DIAZ_full=newConfig('2020.6.23_A__CaMKIIKD_6hbsln_4hgap(1-2)_5uM_48hDIAZ_SPIKES.mat',6);
        sessionConfigs.CaMKIIabKD_keta=newConfig('2020.6.23_B_CaMKIIKD_6hbsln_50uMketa(48h)_SPIKES.mat',6);
        sessionConfigs.Cpt1cKO_Bac_4=newConfig('202.7.9_B_Cpt1c ko_6H BSLN_bac_SPIKES.mat',6);
        sessionConfigs.Stability=newConfig('STABILITY_SPIKES.mat',4);
        sessionConfigs.WT_Bac_0720=newConfig('2020.7.9_A_WT_6hbsln_10uMbac_SPIKES.mat',6);
        sessionConfigs.mini_bac=newConfig('2020-07-14T14-34-25mini_bac_spikes.mat',20);
        sessionConfigs.old_bac=newConfig('2020-07-14T14-47-10Old system_bac.mat',20);
        sessionConfigs.MEAincubator=newConfig('2020-08-02T15-00-22Spikes.mat',20);
        sessionConfigs.CaMKIIabKD_keta2=newConfig('2020.8.12_CamKIIKD_1.2mM_6hbsld_2dKet_washout_spikes0044.mat',6);
        sessionConfigs.Cpt1cKO_Bac5=newConfig('2020.8.16_CPT1cKO_6hbsln_Bac_washout(2h)_spikes.mat',6);
        sessionConfigs.Cpt1aKD_Bac=newConfig('2020-08-20T20-05-27Cpt1aKO_baclofen_6hbsln.mat',6);
        sessionConfigs.Propofol_Bac4=newConfig('2020.8.16_6hbsln_5uM propofol_10uM Bac_spikes.mat',6);
        sessionConfigs.Cpt1cKO_GBZ4=newConfig('2020.8.19_Cpt1cKO_6h bsln_30uM GBZ_spikes.mat',6);
        sessionConfigs.Cpt1cKO_Bac5=newConfig('2020-08-24T15-30-49Cpt1c KO_6hbsln_Bac_spikes.mat',6);
        sessionConfigs.Propofol_newMEA=newConfig('2020-08-21T15-52-10WT_propofol_spikes.mat',6);
        sessionConfigs.Cpt1cKO_Bac5b=newConfig('2020-08-26T16-32-30_Cpt1cKO_all baseline_bac.mat',12);
        sessionConfigs.CaMKIIabKD_GBZ=newConfig('2020-08-30T15-47-54shCaMKII_GBZ_spikes.mat',6);
        sessionConfigs.Paxhiline_Diaz=newConfig('2020-09-03T14-48-27_bsln(11_2h)_30hPaxhiline(2h)_Diaz.mat',6);
        sessionConfigs.CamKIIab_Bac_3=newConfig('2020-09-13T13-14-23MEA1_shCamKII_spikes_6hbsln_Bac.mat',6);
        sessionConfigs.eEF2KKO_Bac4=newConfig('2020-09-13T16-02-32MEA3_eEF2K KO_spikes_6hbsln_Bac.mat',6);
        sessionConfigs.eEF2KKO_10uMGBZ=newConfig('2020-09-14T10-25-20_eEF2KKO_6h bsln_10uM GBZ_NEW_.mat',6);
        sessionConfigs.DPCPX_KeTA_OLD=newConfig('2020.9.14_6h bsln_100nM DPCPX_KETAMINE_spikes.mat',6);
        sessionConfigs.CamKIIab_GBZ_cont=newConfig('2020-09-14T10-29-06_shCamKII_10uM GBZ_contaminated.mat',6);
        sessionConfigs.DPCPX_KeTA_NEW=newConfig('2020-09-21T17-14-21_WT_100nM DPCPX_KETA.mat',6);
        sessionConfigs.CaM1234_keta=newConfig('2020.9.23_CaM1234_6hbsln_50uM keta_spikes.mat',6);
        sessionConfigs.Torin_keta=newConfig('2020.9.23_6HBSLN_100nM Torin(24h)_50uM Keta_spikes.mat',6);
        sessionConfigs.CaM1234_bac=newConfig('2020-09-29T15-46-57_CaM1234_bac.mat',6);
        sessionConfigs.Prozac_test=newConfig('2020.9.29_4h bsln_Prozac.mat',4);
        sessionConfigs.shCamKIIab_keta3=newConfig('2020-10-10T15-14-58_shCamKIIab_6hbsln_keta48h_washout_spikes.mat',6);
        sessionConfigs.shCamKIIa_bac=newConfig('2020.10.6_shCamKIIalpha_6hblsn_bac_spikes.mat',6);
        sessionConfigs.shCamKIIa_keta_bac=newConfig('2020-10-11T21-47-57_shCamKIIalpha_6hbsln_48hKeta_Bac_spikes.mat',6);
        sessionConfigs.shCamKIIb_keta_bac=newConfig('2020-10-11T21-53-37_shCamKIIbeta_6h bsln_keta_Bac_spikes.mat',6);
        sessionConfigs.CPT1cKO_GBZ5=newConfig('2020.10.12_CPT1cKO_6h bslm_10uM GBZ_spikes.mat',6);
        sessionConfigs.CPT1cKO_Bac_cont=newConfig('2020.10.12_CPT1cKO_6h BSL_10uM bac_spike_cont.mat',6);
        sessionConfigs.CPT1cKO_Bac7=newConfig('2020-10-16T17-12-43_CPT1cKO_6h bsln_Bac_SPIKES.mat',6);
        sessionConfigs.Prozac_1to10uM=newConfig('2020.10.19_6hbsln_1uMprozac24h__10uM prozac_10uM bac_spikes.mat',6);
        sessionConfigs.shCamKIIb_keta2=newConfig('2020.10.21_6h_bsln_shCamKIIbeta_50uM Keta_spikes.mat',6);
        sessionConfigs.shCamKIIb_GBZ=newConfig('2020-10-25_shCamKIIbeta_6hbsln_GBZ_spikes.mat',6);
        sessionConfigs.titrKeta_GBZ=newConfig('2020-10-26T13-49-30_0.5uMketa_20uM_50uM_GBZ_spikes.mat',3);
        sessionConfigs.shCamKIIb_keta3=newConfig('2020-10-26T13-52-05_shCamKIIbeta_keta50uM_spikes.mat',6);
        sessionConfigs.shCamKIIb_keta3new=newConfig('2020-10-28T18-57-24_new_shCamKIIbeta_keta50uM_spikes.mat',6);
        sessionConfigs.Prozac_5uM=newConfig('2020-11-08T15-51-28_6hbsln_5uM prozac_bac_spikes.mat',6);
        sessionConfigs.titrKeta_GBZ2=newConfig('2020-11-08T17-35-33_6hrbsln_keta_gbz_spikes.mat',6);
        sessionConfigs.APV50uM_Bac=newConfig('2020-11-09T12-07-23_50uM APV_bac.mat',6);
        sessionConfigs.Torin_Bac2=newConfig('2020-11-15T15-23-19_5x2h bsln_100nM Torin_bac.mat',5);
        sessionConfigs.TeriKeta2=newConfig('2020-11-21T15-17-30_teri_keta.mat',5);
        sessionConfigs.testKD=newConfig('2020-11-25T13-01-40test_KD_24h.mat',5);
        sessionConfigs.testScrmbl=newConfig('2020-11-25T13-05-24Scramble test.mat',5);
        sessionConfigs.CamKIIbetaKD_keta4=newConfig('2020.11.18_Anto_shCamKIIbeta_ketamine_spikes0001.mat',6);
        sessionConfigs.testScrmbl48=newConfig('2020-11-26T15-41-31Scramble test48h.mat',5);
        sessionConfigs.testKD48h=newConfig('2020-11-26T16-38-5448H KD.mat',5);
        sessionConfigs.Scrmbl=newConfig('2020-11-27T17-41-01Scramble_full experiment.mat',5);
        sessionConfigs.CamKIIbKD_KN93=newConfig('2020-11-27T17-45-54_CamKIIbetaKD_full experiment.mat',5);
        sessionConfigs.APV_50uM=newConfig('2020-11-29T19-21-52_50uM APV_SPIKES.mat',5);
        sessionConfigs.Keta_titr_3=newConfig('2020-11-30T12-38-10_ket titration.mat',5);
        sessionConfigs.CamKIIKD_old_system=newConfig('2020.11.18_shCamKIIbeta_72h keta_SPIKES.mat',5);
        sessionConfigs.Prozac_5uM_2=newConfig('2020.11.26_6hbsln_48h_5uM Fluoxetine_bac_SPIKES.mat',6);
        sessionConfigs.CamKIIbetaKD_Keta20uM=newConfig('2020-12-10T16-58-20_CamKIIbetaKD_20uM K_spikes.mat',12);
        sessionConfigs.CamKIIbetaKD_Keta50uM_6=newConfig('2020-12-10T16-50-17_CaMKIIKD_50uM K_spikes.mat',12);
        sessionConfigs.Keta20uM_Bac=newConfig('2020.12.6_Anto_WT20uM K_bac_spikes0001.mat',6);
        sessionConfigs.Ket_titr_1GBZ=newConfig('2020-12-10T16-45-29_K_titration_1uM GBZ_spikes.mat',6);
        sessionConfigs.DAPV_Bac=newConfig('2020-12-10T16-54-33_AP5_Bac.mat',6);
        sessionConfigs.CamKIIbetaKD_Keta10uM_test=newConfig('2020-12-17T10-40-01Test_4h bsln_4h 10uM_spikes.mat',2);
        sessionConfigs.CamKIIalphaKD_10uMK_test=newConfig('2020.12.15_Anto_B_CKIIalphaKD_10uM K_spikes_test.mat',4);
        sessionConfigs.CamKIIbetaKD_Keta10uM_test24h=newConfig('2020-12-17T11-59-22_4h bsln_10uM K_24h_test.mat',4);
        sessionConfigs.DAPV_TEST=newConfig('2020-12-20T14-55-56_APV_test.mat',4);
        sessionConfigs.CamKIIalphaKD_10uMK=newConfig('2020.12.17_12h(1h)bsln_CKIIalphaKD_10uM K_spikes.mat',12);
        sessionConfigs.CamKIIbetaKD_Keta10uM_Bac=newConfig('2020-12-20T14-54-08shCaMKIIbetaKD_12h bsln_10uM K_bac.mat',12);
        sessionConfigs.DAPV_bac2=newConfig('2020-12-27T14-01-23_6hbsln_50uMDAPV_bac_spikes.mat',6);
        sessionConfigs.DAPV_bac3=newConfig('2020.12.20_6hBsln_50uM APV_Bac10M_spikes.mat',6);
        sessionConfigs.MK801_bac=newConfig('2020.12.19_6hrBsln_25uM mk801_Bac10uM_spikes.mat',6);
        sessionConfigs.keta20_Bac_2=newConfig('2020-12-27T16-17-55_12hbsln_20uMKeta_Bac_spikes.mat',12);
        sessionConfigs.MK801_bac2=newConfig('2020.12.29_6h bsln_25uMMMK801_Bac_spikes.mat',6);
        sessionConfigs.keta20_washout=newConfig('2021-01-03T15-05-52_4bsln_20uMKeta_Washout_spikes.mat',8);
        sessionConfigs.MK_Keta=newConfig('2021-01-03T14-55-12_MK25uM_20uMKeta_spikes.mat',12);
        sessionConfigs.APVtitr_keta=newConfig('2021-01-03T15-22-166hBsln_1uM_26h_20uM_28uM_50_24h_20uMKeta_spikes.mat',12);
        sessionConfigs.mida_test=newConfig('2021-01-14T15-14-31mida test1uM.mat',3);
        sessionConfigs.WT_iso_media=newConfig('2021.1.12_WT_6hbsln_iso_spikes.mat',6);
        sessionConfigs.AD_iso_media=newConfig('2021.1.12_AD_6hbsln_iso_spikes.mat',6);
        sessionConfigs.Memantine_Bac1=newConfig('2021-01-17T18-55-06SPIKES_memantine_bac.mat',3);
        sessionConfigs.Memantine_Bac2=newConfig('2021-01-18T11-53-53_MEMANTINE_BAC_SPIKES.mat',3);
        sessionConfigs.Midazolam1uM_Bac=newConfig('2021-01-18T12-03-22_MIDAZOLAM_BAC_SPIKES.mat',3);
        sessionConfigs.DAPV50uM_Ket20uM=newConfig('2021-01-18T11-29-1648H APV_20uM KETA_SPIKES.mat',6);
        sessionConfigs.WT_iso_test2=newConfig('2021.1.17_WT_6hbsln_4%ISO_250MLXMIN_0.1%_spikes.mat',6);
        sessionConfigs.AD_iso_test2=newConfig('2021.1.17_AD_4%ISO_250MLxMIN0.1%_spikes.mat',6);
        sessionConfigs.WT_iso_test3=newConfig('2021.1.20_1DROP ISO_spikes.mat',6);
        sessionConfigs.WT_iso_test2full=newConfig('2021.1.17_WT_6hbsln_4%ISO_250MLXMIN_0.1%_spikes_full.mat',6);
        sessionConfigs.AD_iso_test2full=newConfig('2021.1.17_AD_6hbsln_4%ISO_250MLXMIN_0.1%_spikes_full.mat',6);
        sessionConfigs.Mementine_Bac3=newConfig('2021-02-03T16-32-25_MEMENTINE_BAC_SPIKES.mat',6);
        sessionConfigs.APV_BAC_5=newConfig('2021-02-03T14-57-04_APV_BAC_SPIKES_3.mat',6);
        sessionConfigs.Keta20_GBZ1uM=newConfig('2021-02-03T14-23-39_20uM KETA_1uM GBZ_SPIKES_2.mat',6);
        sessionConfigs.Keta20_GBZ1uM_b=newConfig('2021-02-04T12-17-53_mea4_20uM K_gbz_1uM_spikes.mat',6);
        sessionConfigs.Memantine_GBZ=newConfig('2021-02-09T19-18-52_mea2_MEM30uM_BGZ_1to10_SPIKES.mat',6);
        sessionConfigs.apv_gbz=newConfig('2021-02-09T19-15-20_APV_GBZ1_10_SPIKES_MEA1.mat',6);
        sessionConfigs.Memantine_GBZ_2=newConfig('2021-02-09T19-21-47_MEA4_MEMANT30uM_1to10GBZ_SPIKES.mat',6);
        sessionConfigs.Iso5uL=newConfig('2021.2.7_6hbl_22h5ulISO.mat',6);
        sessionConfigs.Iso7_5uL=newConfig('2021.2.5_6hbl_38h7.5ulISO.mat',6);
        sessionConfigs.shMICU3_Ket=newConfig('2021-02-14T10-23-31_MEA3_shMICU3_SPIKES.mat',6);
        sessionConfigs.CaMKIIaT305D_ket=newConfig('2021-02-14T10-29-29_MEA4_CaMKIIa305D_SPIKES.mat',5);
        sessionConfigs.shMICU3_Ket2=newConfig('2021-02-14T10-09-53MEA1__shMICU3_20keta_spikes.mat',6);
        sessionConfigs.CaMKIIaT305D_ket2=newConfig('2021-02-14T10-16-59_MEA2_CaMKIIaT305D_SPIKES.mat',6);
        sessionConfigs.MCUKO_keta=newConfig('2021-02-22T13-11-18_MCUKO_20uM Ket_spikes.mat',6);
        sessionConfigs.CaMKIIaT305D_ket3=newConfig('2021-02-28T13-53-04_CamKIIaT306D_MEA2_20uM Keta_Spikes.mat',6);
        sessionConfigs.CaMKIIaT305D_ket4=newConfig('2021-02-28T13-53-04_CamKIIaT306D_MEA4_20uM Keta_Spikes.mat',6);
        sessionConfigs.DAPV_GBZ_2=newConfig('2021-03-08T16-29-37_50uMDAPV_10uM GBZ_spikes.mat',6);
        sessionConfigs.DAPV_bac_6=newConfig('2021-03-08T16-33-24_APV_Bac_spikes.mat',6);
        sessionConfigs.Teri100uM_20Keta2=newConfig('2021-03-08T16-36-14_100uM Teri_ket.mat',6);
        sessionConfigs.Keta20uM_Teri100uM=newConfig('2021-03-08T16-44-09_20uM K_100uM Teri.mat',6);
        sessionConfigs.iNeu_H1_bac=newConfig('2021-03-13T12-59-16iNEUR_H1_6H BSLN_10uM BAC_SPIKES.mat',6);
        sessionConfigs.shMICU3_Ket3=newConfig('2021-03-13T12-54-10_shMICU3_20uM KETA_SPIKES.mat',6);
        sessionConfigs.shMICU3_Ket4=newConfig('2021-03-13T13-03-43_shMICU3_6hBSLN_20uMK_SPIKES.mat',6);
        sessionConfigs.MCUff_Keat_test=newConfig('ff spikes.mat',3);
        sessionConfigs.MCkof_Keat_test=newConfig('ko spikes.mat',3);
        sessionConfigs.MCU_KO_Mea3=newConfig('2021-03-17T00-13-31MCU KO_MEA3_SPIKES.mat',6);
        sessionConfigs.H1_GBZ=newConfig('2021-03-17T00-10-41_h1_gbz.mat',6);
        sessionConfigs.MCU_KO_Mea4=newConfig('2021-03-21T13-19-56_MCU KO_MEA4.mat',6);
        sessionConfigs.MCUff_MEA1=newConfig('2021-03-17T00-07-26_MCUff_KETA_SPIKES.mat',6);
        sessionConfigs.MCU_MEA_B=newConfig('2021.3.13_A_MCU KO_20uM KETA_spikes.mat',6);
        sessionConfigs.MCUff_MEA_A=newConfig('2021.3.13_A_MCUff_20uM KETA_first part_spikes.mat',4);
        sessionConfigs.MCUff_MEA_AII=newConfig('2021.3.13_A_MCUff_20uM KETA_IIpart_spikes.mat',4);
        sessionConfigs.MCUff_ket_Bac=newConfig('2021-03-25T15-37-11_MCUff_ketamine_10uM GBZ.mat',6);
        sessionConfigs.DAPV_ket2=newConfig('2021-04-10T14-17-53_MEA1_6h bsln_APV_Ket20uM_spikes.mat',6);
        sessionConfigs.Teri50_ket3=newConfig('2021-04-10T14-23-26_6h bsln_50uM Teri_20uM K_spikes.mat',6);
        sessionConfigs.Teri50_ket4=newConfig('2021-04-10T16-46-02MEA4_5H BSLN_50uM TERI_20uM KET_SPIKES.mat',5);
        sessionConfigs.Mem_Bac4=newConfig('2021-04-11T14-57-31_5H BSLN_30uM MEM_BAC.mat',5);
        sessionConfigs.eEf2KKO_TBOA=newConfig('2021.4.8_EF2KKO_10uM TBOA_SPIKES0001.mat',6);
        sessionConfigs.Mem_Bac5=newConfig('2021.4.8_6h bsln_30uM MEM_BAC_spikes.mat',6);
        sessionConfigs.eEF2KKO_20K_Bac1=newConfig('2021-04-14T20-08-11_eEF2KKO_6h bsln_20uM K_bac_spikes.mat',6);
        sessionConfigs.Teri100_20K_5=newConfig('2021-04-17T18-33-45_6h bsln_100uM Teri_20uM K_spikes.mat',6);
        sessionConfigs.eEF2KKO_20K_Bac2=newConfig('2021-04-17T18-27-29_6hbsln_20uM K_bac_spikes.mat',6);
        sessionConfigs.MK801_bac3=newConfig('2021.4.12_6h bsln_MK-801_10uM BAC_spikes.mat',6);
        sessionConfigs.APV_Ket=newConfig('2021.4.13_6h bsln_48h_50uM APV_50uM TERI_spikes.mat',6);
        sessionConfigs.APV_Ket_3=newConfig('2021-04-25T16-49-37_MEA1_6hbsln_46hDAPV_20uMKeta_spikes.mat',3);
        sessionConfigs.APV_Ket_4=newConfig('2021-04-25T16-55-58_MEA3_6hBSLN_46hDAPV_20uMKET_SPIKES.mat',3);
        sessionConfigs.Mem50_Bac1=newConfig('2021-04-27T12-48-45_6hBsln_50uM MEM_bac_spikes.mat',6);
        sessionConfigs.Mem50_Bac2=newConfig('2021-04-27T12-52-45_6hbsln_50uM MEM_Bac_mea3_spikes.mat',6);
        sessionConfigs.keta20_wash=newConfig('2021-04-28T19-06-17_8h bsln_6h(2h)Ket20uM_Washout.mat',8);
        sessionConfigs.Cnt_5uLISO_Bac=newConfig('2021-05-02T16-26-22_10H BSLN_5uL iso(30_24H)_BAC_SPIKES.mat',10);
        sessionConfigs.APPPS1_5uLISO_Bac=newConfig('2021-05-02T16-59-12_10hBSLN_5uLISO(24h30)_10uM BAC_SPIKES_MEA3.mat',10);
        sessionConfigs.APPPS1_5uLISO_Bac2=newConfig('2021-05-02T16-53-55_APPPS1_12hBSLN_ISO5uL(30_24h)_BAC_SPIKES_MEA1.mat',12);
        sessionConfigs.Cnt_5uLISO_Bac2=newConfig('2021-05-09T16-46-51_WT_7h bsln_12h(30m)_1h_Bac_spikes1.mat',14);
        sessionConfigs.Cnt_5uLISO_Bac2b=newConfig('2021-05-09T16-46-51_WT_7h bsln_12h(30m)_1h_Bac_spikes1b.mat',14);
        sessionConfigs.APPPS1_5uLISO_Bac3=newConfig('2021-05-10T12-05-36_APPPS1_7hbsln_5uLISO(12H_30M_1H)_BAC_SPIKES_2.mat',14);
        sessionConfigs.Cnt_5uLISO_Bac3=newConfig('2021-05-11T11-25-03_WT_7hBSLB_12H_30m_1H 5uL ISO_BAC_SPIKES_4.mat',14);
        sessionConfigs.APPPS1_5uLISO_Bac4=newConfig('2021-05-11T11-17-48_APPPS1_8hBSLN_5uLISO(12_30m_1H)BAC_3_SPIKES.mat',16);
        sessionConfigs.Cntr_5uLISO_Bac4=newConfig('2021-05-16T17-40-13_12hBsln_5uLISOdirect_spikes.mat',12);
        sessionConfigs.Cntr_5uLISO_Bac5=newConfig('2021-05-16T17-43-06_12hBsln_5uLISOinMedia_spikes.mat',12);
        sessionConfigs.Ket1_20_bac=newConfig('2021-05-16T17-46-46_10hBsln_1uMK(5h_1h)_20uMK(1h)_bac(2h)_spikes.mat',10);
        sessionConfigs.Ket20uM_Bac3=newConfig('2021-05-16T17-34-298hBsln_20uMK_Bac_spikes.mat',8);
        sessionConfigs.eEF2KKO_20uMK_3=newConfig('2021-05-22T13-46-20_eEF2KKO_10hBsln_20uMK_spikes.mat',10);
        sessionConfigs.CNTRL_20uMK=newConfig('2021-05-22T13-42-21_Cntrl_10hBSLN_20uMK_spikes.mat',10);
        sessionConfigs.ISO5uLINMEDIA2_BAC=newConfig('2021-05-21T20-35-53_WT_16hbsln_ISO(5uLin50)_14h(1h)_Bac_spikes.mat',16);
        sessionConfigs.cntrl_20uMK_bac3=newConfig('2021-05-24T18-10-21_WT_12hBSLN_20uMK_BAC_SPIKES.mat',12);
        sessionConfigs.eEF2KKO_20uMK_bac3=newConfig('2021-05-24T18-05-32_eEF2KKO_12hBSLN_20uMK_BAC_SPIKES.mat',12);
        sessionConfigs.cntrl_100uMTeri_bac=newConfig('2021-05-25T11-43-23_wt_10h bsln_100uM teri_bac_old_spikes.mat',10);
        sessionConfigs.eEF2KKO_100uMTeri_bac=newConfig('2021-05-23T15-11-03eEF2KKO_12hBSLN_100uMTERI_BAC_SPIKES.mat',12);
        sessionConfigs.cntrl_20uMK_bac4=newConfig('2021-05-30T14-50-09_Cntrl_12hBsln_20uMKeta_Bac_spikes.mat',10);
        sessionConfigs.eEF2KKO_20uMK_4=newConfig('2021-05-30T14-53-49_eEF2KKO_12hBsln_20uMK_Bac_spikes.mat',12);
        sessionConfigs.APV_TERI_2=newConfig('2021-06-06T14-55-54_MEA2_12hBSLN_50uMDAPV_TERI_SPIKES.mat',12);
        sessionConfigs.KET20uM_WASH_2=newConfig('2021-06-06T14-52-32_MEA1_12hBSLN_4h20uM K_WASHOUT(1H)_spikes.mat',12);
        sessionConfigs.KET20uM_WASH_3=newConfig('2021-06-06T15-02-04_MEA4_12hBSLN_4h20uM K_WASHOUT(1H)_spikes.mat',12);
        sessionConfigs.KET20uM_WASH_4=newConfig('2021-06-06T14-58-37_12hBSL_4hKETA_WASHOUT(1h)_SPIKES.mat',12);
        sessionConfigs.ISO5xcent250Mlxmin_A=newConfig('2021.6.6_MEA_A_6h bsln_5%iso_250mLXmin_washout_SPIKES.mat',6);
        sessionConfigs.ISO5xcent250Mlxmin_B=newConfig('2021.6.6_MEA_B_6h bsln_5%iso_250mLXmin_washout_SPIKES.mat',6);
        sessionConfigs.ENBAtestA=newConfig('2021.6.6_MEA_A_4h BSLN_100nMENBA_spikes.mat',4);
        sessionConfigs.ENBAtestB=newConfig('2021.6.6_MEA_B_4h BSL_100nMENBA_spikes0004.mat',4);
        sessionConfigs.ISO_01_1_5_meaA=newConfig('2021-06-10T12-02-01_MEA_A_10hBsln_12h(1h)_(2h)_0.1%ISO_4h1%_1h5%_washout_spikes.mat',10);
        sessionConfigs.ISO_01_1_5_B=newConfig('2021-06-10T12-06-33_MEAB_12hBSLN_0.1_1_5_WASHOUT_SPIKES.mat',12);
        sessionConfigs.Cntr_Teri100uM_2=newConfig('2021-06-13T09-20-00_WT_12hBsln_100uMTERI_Bac_MEA2_spikes.mat',12);
        sessionConfigs.eEF2KKO_Teri100uM_2=newConfig('2021-06-13T09-16-37_eEF2KKO_12hBsln_Teri100uM_Bac_MEA1_spikes.mat',12);
        sessionConfigs.Cntr_Teri100uM_3=newConfig('2021-06-13T10-38-58_WT_6hBSLN_100uM Teri_BAC_MEA4_SPIKES.mat',12);
        sessionConfigs.eEF2KKO_Teri100uM_3=newConfig('2021-06-13T10-34-15_eEF2KKO_6hBsln_100uMTeri_Bac_MEA3_spikes.mat',12);
        sessionConfigs.Cntr_Teri100uM_4=newConfig('2021-06-17T11-47-50_12hBsln_100uMTeri_1h_12h_2h_Bac_spikes.mat',12);
        sessionConfigs.eEF2KKO_Teri100uM_4=newConfig('2021-06-17T12-00-45_eEF2KKO_12hBsln_100uMTeri_12h(1h)_2h_Bac_spikes.mat',12);
        sessionConfigs.eEF2KKO_20uMK_5=newConfig('2021-06-17T11-56-39_eEF2KKO_12hBsln_20uMketa_12(1h)_2h_Bac_spikes.mat',12);
        sessionConfigs.eEF2KKO_20uMK_5c=newConfig('2021-06-17T11-56-39_eEF2KKO_12hBsln_20uMketa_12(1h)_2h_Bac_spikes2.mat',14);
        sessionConfigs.Teri50uM_Test1=newConfig('2021-06-20T11-51-23_6h bsln_50um teri_spikes.mat',6);
        sessionConfigs.Teri50uM_Test2=newConfig('2021-06-20T11-55-07_mea3_test.mat',6);
        sessionConfigs.APV_Teri50uM=newConfig('2021-06-20T18-06-34_8hBSLN_50uMAPV_50uMTERI_SPIKES.mat',8);
        sessionConfigs.ISO_6_5uL_TBOA=newConfig('2021.6.16_B_6.5uLISO_10uM TBOA_SPIKES.mat',6);
        sessionConfigs.APV_Teri50uM_3=newConfig('2021-06-22T14-45-15_MEA4_12hBSLN_APV_TERI50uM_SPIKES.mat',12);
        sessionConfigs.Teri50uM_keta20_1=newConfig('2021-06-22T14-41-11_MEA3_12hBsln_50uMTeri_ket_spikes.mat',12);
        sessionConfigs.Teri100uM_keta20_6=newConfig('2021-06-24T13-03-12_6h bsln_100um teri_ketaspikes.mat',12);
        sessionConfigs.testISO3h=newConfig('3hbsln_3hISO_spikes.mat',3);
        sessionConfigs.testMEA_A=newConfig('6hBsln_24h ISO_spikes.mat',6);
        sessionConfigs.testMEA_B=newConfig('8hbsln_3x2hI_ISO_6hI_washout_SPIKES.mat',8);
        sessionConfigs.MEA_B=newConfig('2021.7.5_8hBsl_3X2hI_ISO_6I_WASHOUT_MEA_B.mat',8);
        sessionConfigs.testTBOA_A=newConfig('ISO2x2h_TBOA(1h)_SPIKES.mat',4);
        sessionConfigs.MEA_A_wt1=newConfig('6hBSL_28hISO_40MIN_TBOA_13x20MIN_1_1h_MERGEDspikes.mat',6);
        sessionConfigs.MEA_A_APP1=newConfig('2021.7.7_MEA_A_APP_6HBsln_6.5uL ISO_spikes.mat',6);
        sessionConfigs.MEA_B_WT2=newConfig('2021.7.7_MEA_B_WT_6HBsln_6.5uL ISO_spike.mat',6);
        sessionConfigs.MEA2_WT3TBOA_1=newConfig('2021-07-09T14-25-27MEA2_WT_3h(1h)Bsl_TBOA(30min)_spikes.mat',3);
        sessionConfigs.MEA4_APP2TBOA_1=newConfig('2021-07-09T14-28-35_MEA4_APP_3hBSLN(1h)_10uMTBOA(30min)_SPIKES.mat',3);
        sessionConfigs.MEA_B_WT2_TBOA=newConfig('ISO2x2h_TBOA(1h)_SPIKES.mat',2);
        sessionConfigs.DHK100uM=newConfig('2021-07-09T19-21-52_DHK 100uM_SPIKES.mat',3);
        sessionConfigs.DHK50uM=newConfig('2021-07-09T19-19-40_50uM DHK_SPIKES.mat',3);
        sessionConfigs.MEA_A_wt1_FULL=newConfig('2021.7.5_A_WT_MERGED EXPERIMENT_spikes.mat',6);
        sessionConfigs.MEA_A_TEST2=newConfig('2hBSLN_ISO_SPIKES.mat',2);
        sessionConfigs.MEA_B_WT3=newConfig('2021.7.7_MEA_B_WT_6HBsln_6.5uL ISO_tboa_full_SPIKES.mat',6);
        sessionConfigs.MEA4_WT_TBOA_2=newConfig('2021-07-12T17-11-15_10.07.21_MEA4_WT_12hBSLN_TBOA_SPIKES.mat',12);
        sessionConfigs.MEA2_APPPS1_TBOA_2=newConfig('2021-07-12T17-07-57_10.07.21_MEA2_APPPS1_12hBSLN_TBOA(30min)_SPIKES.mat',12);
        sessionConfigs.MEA_A_APPPS1_FINAL=newConfig('2021.7.7_MEA_A_APP_6HBsln_11_GAP6h_3X6.5uLISO_tboa_SPIKES.mat',6);
        sessionConfigs.MEA_B_WT_FINAL3=newConfig('2021.7.7_MEA_B_WT_6HBsln_15X6.5uL ISO_tboa_spikes_realone.mat',6);
        sessionConfigs.MEA_4_APPPS1_TBOA_FINAL=newConfig('2021-07-13T12-07-17_7.7.21_MEA4_APP_3hBSLN(1h)_10uMTBOA(1H)_SPIKES.mat',3);
         sessionConfigs.MEA_2_WT_TBOA_FINAL=newConfig('2021-07-13T_7.7.21_MEA2_WT_3h(1h)Bsl_TBOA(1h)_spikes.mat',3);
         sessionConfigs.MEA3_WTB_DHK_ISO=newConfig('2021-07-12T17-56-3512hBSL_DHK(6h_1h_9_2H)100uM_ISO_SPIKES.mat',12);
         sessionConfigs.MEA_A_WT_9_7_4=newConfig('2021.7.9_MEA_A_BALB_WT_24hISO_TBOA_SPIKES.mat',6);
         sessionConfigs.MEA_2_DHK50uM_iso=newConfig('2021-07-13T18-56-31_8.7.21_MEA1_50uMDHK_ISO_SPIKES.mat',12);
         sessionConfigs.MEA_B_APPPS1_10_7_ISO_TBOA=newConfig('2021.7.10_MEA_B_APP_PS1_6hBSL_13hISO_TBOA_SPIKES.mat',6);
         sessionConfigs.MEA2_APPPS1_100uMDHK=newConfig('2021-07-13T12-37-00_12.7.21_APPPS1_4_MEA2_12hBSLN_100uM DHK_SPIKES.mat',12);
         sessionConfigs.MEA_A_WT_11721_ISO_DHK=newConfig('2021.7.11_MEA_A_WT_6hBSLN_13xISO_12x20MIN_1hDHK100mM_SPIKES.mat',6);
         sessionConfigs.MEA3_WT_B_ISO_DHK=newConfig('2021-07-14T15-02-39_MEA3_B_WT_12hBSLN_ISO12_100mM DHK_SPIKES.mat',12);
         sessionConfigs.MEA1_APPPS1_A_ISO_DHK=newConfig('2021-07-14T14-58-39_MEA1_A_appps1_12hBSL_ISO_12_DHK100mM_SPIKES.mat',12);
         sessionConfigs.MEA2_APPPS1_ISO_DHK_2=newConfig('2021-08-01T14-58-18_12.7_MEA2_APPPS1_ISO_DHK_SPIKES.mat',6);
         sessionConfigs.MEA3_WT_ISO_DHK_3=newConfig('2021-08-01T19-11-29_MEA3_WT_ISO_DHK_SPIKES.mat',6);
         sessionConfigs.MEA4_WT_ISO_DHK_4=newConfig('2021-08-01T15-05-46_12.7.21_MEA4_WT_ISO_DHK_SPIKES.mat',6);
         sessionConfigs.MEA1_APPPS1_ISO_DHK_3=newConfig('2021-08-02T10-23-14_APPPS1_MEA1_6hBSLN_ISO_DHK_SPIKES.mat',6);
         sessionConfigs.TrkBFc_Bac1=newConfig('2021.8.17_mea4_12hBSLB_TrkbFc_Bac_spikes.mat',12);
         sessionConfigs.TrkBFc_20Ket1=newConfig('2021.8.18_MEA3_12hrBSLN_TrkBfc_20uMnewK_SPIKES.mat',6);
         sessionConfigs.TrkBFc_TBOA1=newConfig('2021-08-23T12-14-05_MEA2_TrKBFc_TBOA.mat',6);
         sessionConfigs.testISO_bub=newConfig('test_spikes.mat',3);
         sessionConfigs.newKet_TBOA1=newConfig('2021-08-24T10-24-43_12hBSLN_20uMnewKET_TBOA10uM_spikes.mat',12);
         sessionConfigs.newKet_TBOA1b=newConfig('2021-08-24T10-24-43_12hBSLN_20uMnewKET_TBOA10uM_spikes2.mat',12);
         sessionConfigs.testISO2_bub=newConfig('3hbsln_3hISO_spikes.mat',3);
         sessionConfigs.testISO24h_bub2=newConfig('test_4hBsln_iso(4h)_iso(6h1h_2h)_spike.mat',6);
         sessionConfigs.testISO_drop=newConfig('drop test_spikes.mat',4);
         sessionConfigs.testWT_prop=newConfig('2021-08-26T15-42-10_WT_prop_test.mat',3);
         sessionConfigs.testAPP_prop=newConfig('2021-08-26T15-45-17_APP_prop.mat',4);
         sessionConfigs.WT_ISO_complete=newConfig('2021.8.21_WT_6hBSL_100uL100xISO(10)_x2(13)_5uLdrop_spikes0001.mat',6);
         sessionConfigs.test24HAPP_prop=newConfig('2021-08-27T10-55-00_APP_prop24h_spikes.mat',6);
         sessionConfigs.test24HWT_prop=newConfig('2021-08-28T13-44-55_WT_24hprop_spikes.mat',6);
         sessionConfigs.wt_iso_24H_TEST=newConfig('2021.8.26_WT_4hBSLN_6.5uL_Iso_APIKES.mat',4);
         sessionConfigs.APP_iso_24H_TEST=newConfig('2021.8.26_APP_4hBSLN_6.5uL_Iso_SPIKES.mat',4);
         sessionConfigs.TrkBFc_20Ket2_bac=newConfig('2021-08-29T15-54-05_12hBSLN_TrkBFc24h_20uM KET(12h_1h)_BAC_SPIKES.mat',12);
         sessionConfigs.APP_iso_complete=newConfig('2021.8.26_APP_4hBSLN_6.5uL_Iso(6x2h_4h_48h)_100uMTERI_sp.mat',4);
         sessionConfigs.WT_ISO_complete2=newConfig('2021.8.26_WT_4hBSLN_6.5uL(2h_6h_4h48H)_new_Iso_spike.mat',4);
         sessionConfigs.APP_Teri1=newConfig('2021.8.26_APP_8hBSLN_6.5uL_Iso_100uMTERI_spikes.mat',8);
         sessionConfigs.WT_ISO_complete3=newConfig('2021.8.26_WT_8hBSLN_6.5uL_new_Iso(2h)_spikes.mat',8);
         sessionConfigs.WT_prop_compl=newConfig('2021-08-30T16-36-51_WTcomplete_12hBLN_48hprop5uM_spikes.mat',12);
         sessionConfigs.APP_prop_compl=newConfig('2021-08-30T16-41-20_APPcomp_12hBSLN_5uMprop48h_spikes.mat',12);
         sessionConfigs.APP_Teri1_comp=newConfig('2021.8.26_APP_8hBSLN_6.5uL_100uMTERI48h_SPIKES0001.mat',8);
         sessionConfigs.APP_Teri2=newConfig('2021.9.1_APPPS1_8hBSLN_100uM TERI_SPIKES0001.mat',8);
         sessionConfigs.APP_Teri3=newConfig('2021-09-05T16-51-27_APPPS1_12hBSLN_100uM TERI_SPIKES.mat',12);
         sessionConfigs.APP_Bac5=newConfig('2021-09-05T16-58-39_APPPS1_12hBSLN_BAC_SPIKES.mat',12);
         sessionConfigs.WT_Bac7=newConfig('2021-09-06T12-35-27_WT_12hBSLN_BAC_SPIKES.mat',12);
         sessionConfigs.WT_ISO_6uL_6=newConfig('2021.9.1_WT_6hBSLN_6.5uL ISO_SPIKES0001.mat',6);
         sessionConfigs.WT_ISO_5uL_7=newConfig('2021.9.5_WT_5hBSL_5uL ISO_spikes.mat',6);
         sessionConfigs.APPPS1_TERI_3B=newConfig('2021-09-09T10-21-39_mea3_APPPS1_12hBSLN_TERI100uM_SPIKES.mat',12);
         sessionConfigs.NewKET_TEST=newConfig('2021-09-09T10-25-34_12hBSL_20uMnewKETA_SPIKES.mat',12);
         sessionConfigs.TRKBFC_KET_3=newConfig('2021-09-09T10-29-02_MEA4_12hBSLN_TRKBFc(3X2H)_20uM_K(12X1h_2h)_SPIKES.mat',12);
         sessionConfigs.Teri100uM_new=newConfig('2021-09-14T16-55-33_mea4_12hBSLN_100uM TERI_spikes.mat',12);
         sessionConfigs.TRKBFC_KET_4=newConfig('2021-09-15T13-08-08_MEA1_12hBSLN_8hTrkBFc(1h)_20uMk(1h)_SPIKES.mat',12);
         sessionConfigs.Cntrl_GBZ_A1=newConfig('2021.9.20_8hBSLN_30uM GBZ_spike.mat',8);
        sessionConfigs.Cntrl_GBZ_A2=newConfig('2021-09-23T17-05-04_mea4_12hBSL_gbz30uM_spikes.mat',12);
        sessionConfigs.Cntrl_GBZ_A3=newConfig('2021-09-23T17-08-38_mea1_12hBsl_30uM GBZ_spikes.mat',12);
         sessionConfigs.Cntrl_GBZ_A2_SS=newConfig('2021-09-23T17-05-04_mea4_12hBSL_gbz30uM_spikes SORTED.mat',12);
         sessionConfigs.MK_Bac_4=newConfig('2021-09-19T17-16-53_MEA3_MK_BAC_SPIKES.mat',12);
         sessionConfigs.Keta_TBOA_4=newConfig('2021-09-19T17-13-03_MEA2_12hBSLN_20uM K_TBOA_SPIKES.mat',12);
         sessionConfigs.APV_BAC_4=newConfig('2021-10-17T16-37-19_MEA3_12hBSLN_APV_BAC.mat',12);
         sessionConfigs.APV_BAC_5=newConfig('2021-10-17T16-29-01_MEA2_12hbsln_APV_BAC.mat',12);
         sessionConfigs.Keta_TBOA_5=newConfig('2021-10-17T16-40-39_12hBSLN_K_TBOA_spikes.mat',12);
         sessionConfigs.MK_Bac_5=newConfig('2021-10-17T16-18-53_12HBsl_MK_Bac_spikes.mat',12);
         sessionConfigs.TrkBFc_keta_5=newConfig('2021-10-27T13-42-44_8H BSLN_6hTrKBFc_20uM K_SPIKES.mat',16);
         sessionConfigs.TrkBFc_keta_6=newConfig('2021-10-31T16-43-52_MEA2_12hBSLN_TRKB6h_KETA(12h_1H)_SPIKES.mat',12);
         sessionConfigs.TrkBFc_keta_7=newConfig('2021-10-31T16-37-45_MEA1_12hBSLN_TRKB6h_KETA(12h_1H)_SPIKES.mat',12);
         sessionConfigs.TrkBFc_keta_8=newConfig('2021-10-31T16-48-05_MEA3_12hBSLN_TRKBandKETA(12h_1H)_SPIKES.mat',12);
         sessionConfigs.TrkBFc_keta_9=newConfig('2021-10-31T16-52-29_MEA4_12hBSLN_TRKBFCandKETA(12h1h)_SPIKES.mat',12);
         sessionConfigs.TalT_1test=newConfig('2021-11-04T11-51-11_mea2_8h bsln_TalT_test_spikes.mat',8);
         sessionConfigs.TalT_2test=newConfig('2021-11-04T11-53-20_mea3_8h Bsln_TalT test.mat',8);
         sessionConfigs.Noradrenalin_1test=newConfig('2021-11-04T11-48-19_mea1_test_8h BSL_norep_spikes.mat',8);
         sessionConfigs.Noradrenalin_2test=newConfig('2021-11-04T11-56-14_mea4_8h bsln_norep test_spikes.mat',8);
         sessionConfigs.TalT_1=newConfig('2021-11-08T13-03-57_mea2_8h bsln_1uMTalT_spikes.mat',8);
         sessionConfigs.TalT_2=newConfig('2021-11-08T13-06-23_mea3_8h Bsln_1uM_TalT_spikes.mat',8);
         sessionConfigs.Noradrenalin_1=newConfig('2021-11-08T13-01-51_mea1_8h BSL_10uM NA_spikes.mat',8);
         sessionConfigs.Noradrenalin_2=newConfig('2021-11-08T13-09-32_mea4_8h bsln_10uM NA_spikes.mat',8);
         sessionConfigs.MCH_1=newConfig('2021-11-21T11-01-46_12hBSL_1uM MHC_spikes.mat',12);
         sessionConfigs.MCH_2_BAC=newConfig('2021-11-21T11-24-21_MEA3_1uM MHC_BAC_SPIKES.mat',12);
         sessionConfigs.TrkBFc_keta_10=newConfig('2021-11-28T11-41-12_MEA1_12hBSL_TRKBFcKETA_SPIKES.mat',12);
          sessionConfigs.TrkBFc_Bac_2=newConfig('2021-11-28T11-47-21_MEA3_8hBSL_18hTrkBFc_TRKBFCBAC_SPIKES.mat',8);
          sessionConfigs.Keta_TBOA_6=newConfig('2021-11-28T11-53-25MEA4_12hBSL_20uM KETA_26H_TBOA_SPIKES.mat',12);
          sessionConfigs.eEF2KO_KETA_6=newConfig('2021-12-03T09-17-06_MEA1_eEF2KKO_KET_SPIKES.mat',12);
          sessionConfigs.eEF2KO_TBOA_3=newConfig('2021.11.30_B_eEF2KKO_6hBSL_10uM TBOA_spikes.mat',6);
         sessionConfigs.eEF2KO_TBOA_4=newConfig('2021.11.30_A_eEF2KKO_6hBSL_10uM TBOA_sp.mat',6);
         sessionConfigs.WT_TBOA_1221=newConfig('2021-12-06T11-02-57_MEA2_12hBSLN_TBOA_SPIKES.mat',12);
         sessionConfigs.TrkBFc_TBOA2=newConfig('2021-12-06T11-06-11_12hBSLN_TRKBFcANDTBOA_SPIKES_one.mat',12);
         sessionConfigs.TrkBFc_TBOA_new=newConfig('2021-12-06T11-09-55_MEA1_6hBSLN_20hTRKBFc_TBOAANDTRKB_SPIKES_new.mat',12);
         sessionConfigs.Teri_test12=newConfig('2021-12-09T12-15-14_test_teri.mat',3);
         sessionConfigs.Trk_test12=newConfig('2021-12-09T12-18-17_test trkb.mat',3);
         sessionConfigs.Teri_cntrl=newConfig('2021-12-14T17-00-54_mea1_12hBsln_50uMTERi_spikes.mat',12);
         sessionConfigs.TrkBTeri_1=newConfig('2021-12-14T17-10-18_mea2_12hBsn_trkbteri(12h1h)_spikes.mat',12);
         sessionConfigs.TrkBTeri_2=newConfig('2021-12-14T17-14-14_mea3_12hBsln_TrkbTeri_spikes.mat',12);
         sessionConfigs.TrkBFc_TBOA3=newConfig('2021.12.15_B_6hBSLN_trkbtboa_spikes0001.mat',6);
         sessionConfigs.TBOA_1512wt=newConfig('2021.12.15_A_problem_6hBSLN_10uM tboa_spikes.mat',6);
         sessionConfigs.Teri_cntrl_3=newConfig('2021-12-18T15-24-16_2_12hBSL_50uMTeri_spikes.mat',12);
         sessionConfigs.TrkBTeri_3=newConfig('2021-12-18T15-21-03_1_12hBSL_TrkBTeri50uM_spikes.mat',12);
         sessionConfigs.TrkB_Bac_3=newConfig('2021-12-18T15-27-30_3_12hBSL_trkb_bac_spikes.mat',12);
         sessionConfigs.Serine10uM_2=newConfig('2021-12-18T15-31-17_4_12hBsl_10uMDserine_spikes.mat',12);
         
         %%spike sorted file for NMDAR project%%
         %CONTROL 20uM KETAMINE%
         sessionConfigs.EXP1_2712_BAC=newConfig('2020-12-27T16-17-55_12hbsln_20uMKeta_Bac_spikes.mat',12);
         sessionConfigs.EXP2_165_1_BAC=newConfig('2021-05-16T17-34-298hBsln_20uMK_Bac_spikes SORTED.mat',16);
         sessionConfigs.EXP3_165_2_BAC=newConfig('2021-05-16T17-46-46_10hBsln_1uMK(5h_1h)_20uMK(1h)_bac(2h)_spikes.mat',10);
         sessionConfigs.EXP4_9921=newConfig('2021-09-09T10-25-34_12hBSL_20uMnewKETA_SPIKES.mat',12);
         sessionConfigs.EXP5_199_TBOA=newConfig('2021-09-19T17-13-03_MEA2_12hBSLN_20uM K_TBOA_SPIKES SORT.mat',12);
         sessionConfigs.EXP6_1719_TBOA=newConfig('2021-10-17T16-40-39_12hBSLN_K_TBOA_spikes.mat',12);
         %eEF2K KO 20uM KETAMINE%
         sessionConfigs.EXP1_14_4_21=newConfig('2021-04-14T20-08-11_eEF2KKO_6h bsln_20uM K_bac_spikes_ SORT.mat',6);
         sessionConfigs.EXP2_17_4_21=newConfig('2021-04-17T18-27-29_6hbsln_20uM K_bac_spikes_sorted.mat',6);
         %sessionConfigs.EXP3_22_5_21=newConfig('2021-05-22T13-46-20_eEF2KKO_10hBsln_20uMK_spikes sorted.mat',10);
         sessionConfigs.EXP4_26_5_21=newConfig('2021-05-24T18-05-32_eEF2KKO_12hBSLN_20uMK_BAC_SPIKES_sort.mat',12);
         sessionConfigs.EXP5_30_5_21=newConfig('2021-05-30T14-53-49_eEF2KKO_12hBsln_20uMK_Bac_spikes_sorted.mat',12);
         
         
         
end

function config = newConfig(filePath, nBaselineHours)
    config = struct;
    config.filePath = filePath;
    config.nBaselineHours = nBaselineHours;
end