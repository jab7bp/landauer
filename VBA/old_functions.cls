Public Function calc_photon_Cu_SDE(PL As Double, Al As Double, Cu As Double)
    Dim PL_Cu As Double
    Dim Al_Cu As Double
    Dim photon_Cu_SDE As Double
    Dim Fact As Double
    
    Dim c0, c1, c2, c3, c4, c5, c6 As Double
    c0 = -6.70409912268555E-02
    c1 = 3.55310017227627E-03
    c2 = -4.13560240497645E-05
    c3 = 0.852762573064795
    c4 = 4.04892595252728E-02
    c5 = 0.014925439998092
    
    PL_Cu = PL / Cu
    Al_Cu = Al / Cu
    
    If (Al_Cu > 15#) And (PL_Cu > 35#) Then
        photon_Cu_SDE = Cu / 0.03
    Else
        Fact = (1# + (c0 * (Al_Cu)) + _
                    (c1 * ((Al_Cu) ^ 2)) + _
                    (c2 * ((Al_Cu) ^ 3)))
        If Fact <> 0# Then
            photon_Cu_SDE = ((c3 + (c4 * (Al_Cu)) + _
                (c5 * ((Al_Cu) ^ 2))) / Fact) * Cu
        Else
        End If
    End If
    
    calc_photon_Cu_SDE = photon_Cu_SDE
End Function

Function calc_radiation_quality(OW As Double, PL As Double, Al As Double, Cu As Double, Optional return_type As Integer = 1, Optional rad_qual_mixed_input As Double = -99)
    'Start by looking for betas and classifying that and then move to photons
    Dim final_rad_quality As Integer
    Dim final_rad_quality_string As String
    'final_rad_quality = 0 --> Photon
    'final_rad_quality = 1 --> Beta
    
    Dim beta_indicator As Integer
    beta_indicator = jboyd_calc_mixed_beta(OW, PL, Al, Cu)
    
    Dim photon_source As Double
    Dim photon_source_string As String
    Dim energy As Double
            
    photon_source = 0
    photon_source_string = ""
    'photon_source = 0 --> No photons
    'photon_source = 1 --> Low CVs, High Energy, PH
    'photon_source = 2 --> Mid-energy, PM
    'photon_source = 3 --> Low-energy, PL
    
    'Beta Low = 0
    'Beta High = 1
    
    'return_type = 1 --> Return an integer
    'return_type = 0 --> Return String
    
    Dim beta_low_high As Integer
    beta_low_high = 0
    Dim beta_low_high_string As String
    beta_low_high_string = ""
    
    'Mixed sources
    Dim mixed_beta As Integer
    'If no input (mixed_beta_input) is provided, then we look up the mixed_beta condition.
    'mixed_beta = 0 --> NOT mixed
    'mixed_beta = 1 --> Mixed
    
    Dim signal_sum As Double
    
    Dim OW_PL As Double
    Dim OW_Al As Double
    Dim OW_Cu As Double
    Dim PL_Cu As Double
    Dim OW_PL_cuCorr As Double
    
    OW_PL = OW / PL
    OW_Al = OW / Al
    OW_Cu = OW / Cu
    PL_Cu = PL / Cu
    
    If (PL = Cu) Then
'        OW_PL_cuCorr = 9999
        PL = PL * 1.00000001
    Else
        OW_PL_cuCorr = (OW - Cu) / (PL - Cu)
    End If
    
    signal_sum = OW + PL + Al + Cu
    
    Dim OW_PL_SumScaled As Double
    Dim OW_Al_SumScaled As Double
    Dim OW_Cu_SumScaled As Double
    Dim PL_Cu_SumScaled As Double
    Dim OW_PL_CuCorr_SumScaled As Double
    
    OW_PL_SumScaled = 1000 * OW_PL / signal_sum
    OW_Al_SumScaled = 1000 * OW_Al / signal_sum
    OW_Cu_SumScaled = 1000 * OW_Cu / signal_sum
    PL_Cu_SumScaled = 1000 * PL_Cu / signal_sum
    OW_PL_CuCorr_SumScaled = 1000 * OW_PL_cuCorr / signal_sum
    
    Dim OW_PL_max As Double
    Dim OW_Al_max As Double
    Dim OW_Cu_max As Double
    Dim PL_Cu_max As Double
    Dim OW_PL_CuCorr_max As Double
    
    Dim tol_OW_PL As Double
    Dim tol_OW_Al As Double
    Dim tol_OW_Cu As Double
    Dim tol_PL_Cu As Double
    Dim tol_OW_PL_CuCorr As Double
    
    tol_OW_PL = 0
    tol_OW_Al = 0
    tol_OW_Cu = 0
    tol_PL_Cu = 0
    tol_OW_PL_CuCorr = 0
    
    'OW_PL_max = (1 + tol_OW_PL) * 7
    OW_Al_max = (1 + tol_OW_Al) * 17
    'OW_Cu_max = (1 + tol_OW_Cu) * 18
    'PL_Cu_max = (1 + tol_PL_Cu) * 15
    'OW_PL_CuCorr_max = (1 + tol_OW_PL_CuCorr) * 7
    
    'If (OW_Al_SumScaled > OW_Al_max) Then
    '    beta_low_high = 1
    '    beta_low_high_string = "BH"
    'Else
    '    beta_low_high = 0
    '    beta_low_high_string = "BL"
    'End If
    
    'Calculate Effective Energy
    If (PL_Cu) <= 1.05 Then
        energy = 662#
    Else
        energy = (516.651714671778 + (-4162.31362737238 * Math.Log(PL_Cu)) + _
            (13675.5094026176 * Math.Log(PL_Cu) * Math.Log(PL_Cu)) + _
            (-544.223168651357 * (Math.Log(PL_Cu)) ^ 3)) / _
            (1 + (7.57260700736831 * Math.Log(PL_Cu)) + _
            (35.3741528505647 * Math.Log(PL_Cu) * Math.Log(PL_Cu)) + _
            (121.965738601505 * (Math.Log(PL_Cu)) ^ 3))
    End If

'BETA
    If (rad_qual_mixed_input < 0) Then
        mixed_beta = jboyd_calc_mixed_beta(OW, PL, Al, Cu)
    Else
        mixed_beta = rad_qual_mixed_input
    End If
    
    If (OW_PL_CuCorr_SumScaled = 0) Then
        beta_low_high = 1
        If (mixed_beta) Then
            beta_low_high_string = "BHM"
        Else
            beta_low_high_string = "BH"
        End If
    ElseIf ((OW_Al_SumScaled / OW_PL_CuCorr_SumScaled > 1.5) And (OW_Al_SumScaled > 0.2) And (OW_Cu > 45)) Or (OW_PL_CuCorr_SumScaled <= 0) Or (OW_PL > 20) Then
        beta_low_high = 0
        If (mixed_beta) Then
            beta_low_high_string = "BLM"
        Else
            beta_low_high_string = "BL"
        End If
    Else
        If (mixed_beta) Then
            beta_low_high = 1
            beta_low_high_string = "BHM"
        Else
            beta_low_high = 77
            beta_low_high_string = "BM"
        End If
    End If

    'If (OW_PL > 1.25) And (OW_Al > 1.25) And (OW_Cu > 1.25) And (PL_Cu < 20) And (OW_PL_cuCorr > 20) Then
    '    beta_low_high = 0
    '    beta_low_high_string = "BL"
    'ElseIf (OW_Al > 1.25) And (OW_Cu > 1.25) And (PL_Cu < 17) And (OW_PL_cuCorr > 1.25) And (OW_PL_cuCorr <= 20) Then
    '    beta_low_high = 1
    '    beta_low_high_string = "BH"
    'Else
    '    beta_low_high = 0
    '    beta_low_high_string = "BL"
    'End If
'------------------------
'Photon
    If (beta_indicator = 0) Then
        If (2 < OW <= 20) And (2 < PL <= 20) And (2 < Cu <= Cu) And (2 < Al <= 20) Then
            photon_source = 1
            photon_source_string = "PH"
            beta_indicator = 0
        ElseIf (energy <= 40#) Then
                photon_source = 3
                photon_source_string = "PL"
        ElseIf (40# < energy <= 200) Then
                photon_source = 2
                photon_source_string = "PM"
        Else
                photon_source = 1
                photon_source_string = "PH"
        End If
    End If
    
'----------------------
    'Cases for final output (Either photon or beta)
    If (photon_source = 0) Then
        final_rad_quality = 0
        final_rad_quality_string = photon_source_string
    Else
        final_rad_quality = 1
        final_rad_quality_string = beta_low_high_string
    End If
    
    If (return_type) Then
        'calc_radiation_quality = beta_low_high
        calc_radiation_quality = final_rad_quality
    Else
        'calc_radiation_quality = beta_low_high_string
        calc_radiation_quality = final_rad_quality_string
    End If

End Function

Function old_calc_beta_SDE(dblOW As Double, dblPL As Double, dblAl As Double, dblCu As Double)
    Dim dblBetaSDE As Double
    'Calculate Beta DDE & Beta SDE
    If dblOWPL > 1.25 And dblOWAl > 1.25 And dblOWCu > 1.25 And dblPLCu < 20 And dblOWCU_PLCU_Ratio > 20 Then
        dblBetaSDE = (2.3 * dblOW) + (-2.3 * dblPL) + (1 * dblCu)
    ElseIf dblOWAl > 1.3 And dblOWAl < 3# And dblOWCU_PLCU_Ratio > 1.25 Then
        dblBetaSDE = (1.4624 * dblOW) + (-1.4983 * dblPL) + (0.3246 * dblAl) + (0.6824 * dblCu)
    Else
        dblBetaSDE = (1.4624 * dblOW) + (-1.4983 * dblPL) + (0.3246 * dblAl) + (0.6824 * dblCu)
    End If
    
    old_calc_beta_SDE = dblBetaSDE
End Function

Function old_calc_beta_DDE(dblOW As Double, dblPL As Double, dblAl As Double, dblCu As Double)
    Dim dblBetaDDE As Double
    'Calculate Beta DDE & Beta SDE
    If dblOWPL > 1.25 And dblOWAl > 1.25 And dblOWCu > 1.25 And dblPLCu < 20 And dblOWCU_PLCU_Ratio > 20 Then
        dblBetaDDE = dblCu
    ElseIf dblOWAl > 1.3 And dblOWAl < 3# And dblOWCU_PLCU_Ratio > 1.25 Then
        dblBetaDDE = (-0.3083 * dblOW) + (0.7495 * dblPL) + (-0.3656 * dblAl) + (0.9991 * dblCu)
    Else
        dblBetaDDE = dblCu
    End If
    
    old_calc_beta_DDE = dblBetaDDE
End Function

Function old_calc_photon_Cu_SDE(dblPL As Double, dblAl As Double, dblCu As Double)
    Dim dblFact, dblCuSDE As Double
    Dim dblAlCu As Double
    dblAlCu = dblAl / dblCu
    
    'Calculate Photon Cu SDE
    If (dblAlCu) > 15# And (dblPL / dblCu) > 35# Then
        dblCuSDE = dblCu / 0.03
    Else
        dblFact = (1# + (-6.70409912268555E-02 * (dblAlCu)) + _
            (3.55310017227627E-03 * ((dblAlCu) ^ 2)) + _
            (-4.13560240497645E-05 * ((dblAlCu) ^ 3)))
        If dblFact <> 0# Then
            dblCuSDE = ((0.852762573064795 + (4.04892595252728E-02 * (dblAlCu)) + _
                (0.014925439998092 * ((dblAlCu) ^ 2))) / dblFact) * dblCu
        Else
        End If
    End If
    old_calc_photon_Cu_SDE = dblCuSDE
End Function

Function old_calc_photon_PL_SDE(dblPL As Double, dblCu As Double)
    Dim dblPLSDE As Double
    'Calculate Photon PL SDE
    dblPLSDE = ((1.03317386542118 + _
        (0.899082803473556 * Math.Log(dblPL / dblCu)) + _
        (4.14096846442675E-02 * (Math.Log(dblPL / dblCu) ^ 2)) + _
        (0.128685795322186 * (Math.Log(dblPL / dblCu) ^ 3)) + _
        (-0.165167404396516 * (Math.Log(dblPL / dblCu) ^ 4)) + _
        (2.85842168797508E-02 * (Math.Log(dblPL / dblCu) ^ 5))) / _
        (1# + (1.25183687443143 * Math.Log(dblPL / dblCu)) + _
        (5.04014189193258 * (Math.Log(dblPL / dblCu) ^ 2)) + _
        (-2.9289952309639 * (Math.Log(dblPL / dblCu) ^ 3)) + _
        (0.396407968752253 * (Math.Log(dblPL / dblCu) ^ 4)) + _
        (3.04114687307896E-04 * (Math.Log(dblPL / dblCu) ^ 5)))) * dblPL
        
    old_calc_photon_PL_SDE = dblPLSDE
End Function

Function old_calc_average_photon_SDE(dblPL As Double, dblAl As Double, dblCu As Double)
    Dim dblCuSDE, dblPLSDE, dblAvgSDE As Double
    dblCuSDE = old_calc_photon_Cu_SDE(dblPL, dblAl, dblCu)
    dblPLSDE = old_calc_photon_PL_SDE(dblPL, dblCu)
    'Calculate Photon Average SDE
    dblAvgSDE = (dblCuSDE + dblPLSDE) / 2
    
    old_calc_average_photon_SDE = dblAvgSDE
End Function

Function old_calc_photon_Cu_DDE(dblAl As Double, dblCu As Double)
    Dim dlbCuDDE, dblAlCu As Double
    
    dblAlCu = dblAl / dblCu
    'Calculate Photon Cu DDE
    dblCuDDE = ((0.729750759967566 + (0.129198337923333 * ((dblAlCu) ^ 2))) ^ (0.5)) * dblCu
    
    old_calc_photon_Cu_DDE = dblCuDDE
End Function

Function old_calc_photon_PL_DDE(dblPL As Double, dblCu As Double)
    Dim dblPLDDe As Double
    'Calculate Photon PL DDE
    dblPLDDe = ((2.27875800359032 + _
        (-0.160519484455387 * (dblPL / dblCu)) + _
        (0.444931311949578 * ((dblPL / dblCu) ^ 2)) + _
        (-1.39776305022687E-02 * ((dblPL / dblCu) ^ 3)) + _
        (1.25778369690611E-04 * ((dblPL / dblCu) ^ 4))) / _
        (1# + (-6.65137516928388E-02 * (dblPL / dblCu)) + _
        (1.49910566096546 * ((dblPL / dblCu) ^ 2)) + _
        (-4.14191308266288E-02 * ((dblPL / dblCu) ^ 3)) + _
        (2.2610007373255E-04 * ((dblPL / dblCu) ^ 4)) + _
        (2.1999287177623E-06 * ((dblPL / dblCu) ^ 5)))) * dblPL
    
    old_calc_photon_PL_DDE = dblPLDDe
End Function

Function old_calc_average_photon_DDE(dblPL As Double, dblAl As Double, dblCu As Double)
    'PREVIOUS ALGORITHM FOR calc_average_photon_DDE
    
    Dim dblAvgDDE, dblCuDDE, dblPLDDe As Double
    
    dblCuDDE = old_calc_photon_Cu_DDE(dblAl, dblCu)
    dblPLDDe = old_calc_photon_PL_DDE(dblPL, dblCu)

    'Calculate Photon Average DDE
    dblAvgDDE = (dblCuDDE + dblPLDDe) / 2
    
    old_calc_average_photon_DDE = dblAvgDDE
End Function

Function old_beta_indicator_calc(OW As Double, PL As Double, Al As Double, Cu As Double)
    'PREVIOUS ALGORITHM FOR BETA_INDICATOR
    
    Dim BetaIndicator As Integer
    Dim OW_PL As Double
    Dim OW_Al As Double
    Dim OW_Cu As Double
    Dim PL_Cu As Double
    Dim OW_PL_cuCorr As Double
    
    OW_PL = OW / PL
    OW_Al = OW / Al
    OW_Cu = OW / Cu
    PL_Cu = PL / Cu
    
    If (PL = Cu) Then
'        OW_PL_cuCorr = 0.000001
        PL = PL * 1.00000001
    Else
        OW_PL_cuCorr = (OW - Cu) / (PL - Cu)
    End If
    'Calculate Beta Indicator
     If OW_PL > 1.25 And _
        OW_Al > 1.25 And _
        OW_Cu > 1.25 And _
        PL_Cu < 20 And _
        OW_PL_cuCorr > 20 Then
        BetaIndicator = 1
     ElseIf OW_Al > 1.25 And _
            OW_Cu > 1.25 And _
            PL_Cu < 17 And _
            OW_PL_cuCorr > 1.25 And _
            OW_PL_cuCorr <= 20 Then
            BetaIndicator = 1
     Else
         BetaIndicator = 0
     End If
     
     old_beta_indicator_calc = BetaIndicator
End Function

Function old_calc_effective_energy(dblPL As Double, dblCu As Double)
    'PREVIOUS ALGORITHM FOR EFFECTIVE ENERGY
    Dim dblEnergy As Double
    Dim dblPLCu As Double

    dblPLCu = dblPL / dblCu

'Calculate Effective Energy
    If dblPLCu <= 1.05 Then
        dblEnergy = 662#
    Else
        dblEnergy = (516.651714671778 + (-4162.31362737238 * Math.Log(dblPL / dblCu)) + _
            (13675.5094026176 * Math.Log(dblPL / dblCu) * Math.Log(dblPL / dblCu)) + _
            (-544.223168651357 * (Math.Log(dblPL / dblCu)) ^ 3)) / _
            (1 + (7.57260700736831 * Math.Log(dblPL / dblCu)) + _
            (35.3741528505647 * Math.Log(dblPL / dblCu) * Math.Log(dblPL / dblCu)) + _
            (121.965738601505 * (Math.Log(dblPL / dblCu)) ^ 3))
    End If
    
    old_calc_effective_energy = dblEnergy

End Function


Function old_radiation_quality(dblOW As Double, dblPL As Double, dblCu As Double, dblAl As Double)
        'PREVIOUS ALGORITHM FOR RADIATION QUALITY

        Dim blnLow As Boolean
        Dim dblRatio2 As Double
        Dim rstrRadQuality, PGC_RAD_QLTY_PH, PGC_RAD_QLTY_PL, PGC_RAD_QLTY_PM, PGC_RAD_QLTY_P As String
        Dim PGC_RAD_QLTY_BL, PGC_RAD_QLTY_BH As String
        Dim dblOWCU_PLCU_Ratio As Double
        Dim dblOWPL, dblOWCu, dblOWAl, dbPLCu As Double
        
        dblOWPL = dblOW / dblPL
        dblOWCu = dblOW / dblCu
        dblOWAl = dblOW / dblAl
        dblPLCu = dblPL / dblCu

        PGC_RAD_QLTY_PH = "PH"
        PGC_RAD_QLTY_PL = "PL"
        PGC_RAD_QLTY_PM = "PM"
        PGC_RAD_QLTY_P = "P"
        
        PGC_RAD_QLTY_BH = "BH"
        PGC_RAD_QLTY_BL = "BL"
        
        Dim dblEnergy As Double
        Dim blnBetaIndicator As Boolean
        blnBetaIndicator = old_beta_indicator_calc(dblOW, dblPL, dblAl, dblCu)
        dblEnergy = old_calc_effective_energy(dblPL, dblCu)
        
        If (dblPL = dblCu) Then
'            dblOWCU_PLCU_Ratio = 9999
            PL = PL * 1.00000001
        Else
            dblOWCU_PLCU_Ratio = (dblOW - dblCu) / (dblPL - dblCu)
        End If
        
        'Calculate Radiation Quality
        If dblOW > 2 And dblOW <= 20 And _
           dblPL > 2 And dblPL <= 20 And _
           dblCu > 2 And dblCu <= 20 And _
           dblAl > 2 And dblAl <= 20 Then
            'low CVs
            rstrRadQuality = PGC_RAD_QLTY_PH
            blnBetaIndicator = False

        ElseIf blnBetaIndicator = False Then
            'photon only
            If dblEnergy <= 40# Then
                rstrRadQuality = PGC_RAD_QLTY_PL
            ElseIf dblEnergy > 40# And dblEnergy <= 200# Then
                rstrRadQuality = PGC_RAD_QLTY_PM
            Else
                rstrRadQuality = PGC_RAD_QLTY_PH
            End If

        Else
            'beta
            If dblOWPL > 1.25 And dblOWAl > 1.25 And dblOWCu > 1.25 And dblPLCu < 20 And dblOWCU_PLCU_Ratio > 20 Then
                rstrRadQuality = PGC_RAD_QLTY_BL
            ElseIf dblOWAl > 1.25 And dblOWCu > 1.25 And dblPLCu < 17 And dblOWCU_PLCU_Ratio > 1.25 And dblOWCU_PLCU_Ratio <= 20 Then
                rstrRadQuality = PGC_RAD_QLTY_BH
            Else
                rstrRadQuality = PGC_RAD_QLTY_BL
            End If
        End If
        old_radiation_quality = rstrRadQuality
End Function

