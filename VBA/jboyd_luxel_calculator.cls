Function Largest(ParamArray A() As Variant) As Variant
'returns the largest element of list
'List is supposed to be consistent: all nummbers or all strings
'e.g:   largest(2,6,-9,7,3)         -> 7
'       largest("d", "z", "c", "x") -> "z"

    Dim result As Variant
    Dim i As Integer

    result = Null

    For i = LBound(A) To UBound(A)
        If result > A(i) Then
            'nothing to do. This construct will properly handle null values
        Else
            result = A(i)
        End If
    Next i
    Largest = result
End Function

Function Smallest(ParamArray A() As Variant) As Variant
'returns the largest element of list
'List is supposed to be consistent: all nummbers or all strings
'e.g:   largest(2,6,-9,7,3)         -> 7
'       largest("d", "z", "c", "x") -> "z"

    Dim result As Variant
    Dim i As Integer

    result = Null

    For i = LBound(A) To UBound(A)
        If result < A(i) Then
            'nothing to do. This construct will properly handle null values
        Else
            result = A(i)
        End If
    Next i
    Smallest = result
End Function

Function SmallestPositive(ParamArray A() As Variant) As Variant
'returns the largest element of list
'List is supposed to be consistent: all nummbers or all strings
'e.g:   largest(2,6,-9,7,3)         -> 7
'       largest("d", "z", "c", "x") -> "z"

    Dim result As Variant
    Dim i As Integer

    result = Null

    For i = LBound(A) To UBound(A)
        If A(i) > 0 Then
            If result < A(i) Then
                'nothing to do. This construct will properly handle null values
            Else
                result = A(i)
            End If
        End If
    Next i
    SmallestPositive = result
End Function

Public Function jboyd_calc_effective_energy(PL As Double, Cu As Double)
    Dim energy As Double
    Dim PL_Cu As Double
    
    PL_Cu = PL / Cu
    
    If (PL_Cu <= 1.05) Then
        energy = 662
    Else
        energy = (516.651714671778 + (-4162.31362737238 * Math.Log(PL_Cu)) + _
                    (13675.5094026176 * Math.Log(PL_Cu) * Math.Log(PL_Cu)) + _
                    (-544.223168651357 * (Math.Log(PL_Cu)) ^ 3)) / _
                    (1 + (7.57260700736831 * Math.Log(PL_Cu)) + _
                    (35.3741528505647 * Math.Log(PL_Cu) * Math.Log(PL_Cu)) + _
                    (121.965738601505 * (Math.Log(PL_Cu)) ^ 3))
    End If
    
    jboyd_calc_effective_energy = energy
    
End Function

Function jboyd_calc_photons_combined_SDE(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c0, c1, c2, c3, c4, Sum, combined_SDE As Double
    
    Sum = OW + PL + Al + Cu
    
    c1 = 0.406
    c2 = 0.764
    c3 = -1.47
    c4 = 1.285

    combined_SDE = c1 * OW + c2 * PL + c3 * Al + c4 * Cu
    
    jboyd_calc_photons_combined_SDE = combined_SDE
End Function

Function jboyd_calc_mixed_combined_SDE(OW As Double, PL As Double, Al As Double, Cu As Double)

    'mixed_beta_indicator = 0 --> NOT MIXED
    'mixed_beta_indicator = 1 --> MIXED
    Dim mixed_beta_indicator As Integer
    Dim c0, c1, c2, c3 As Double
    
    Dim OW_Cu, PL_Cu As Double

    If (PL = Cu) Then
        PL = 1.000001 * PL
    End If
    PL_Cu = PL / Cu
    OW_Cu = OW / Cu
    
    'I think we can get a more accurate calculation if we can distinguish between M30 and the others.
    'Let's try to isolate M30 by checking OW/Cu andPL/Cu
    'We will use linear functions from SDE vs (OW + PL + Al + Cu):
    Dim m30_m, m30_b, Sum As Double
    Dim all_m, all_b, combined_SDE, all_sum As Double
    
'Standard Coefficients Method:

    If (OW_Cu > 22) And (PL_Cu > 22) Then  'M30
        
        'Coefficient Method:
        c0 = 0.43101
        c1 = 0.0619596
        c2 = -0.533692
        c3 = 9.04421

        combined_SDE = 1.08 * (c0 * OW + c1 * PL + c2 * Al + c3 * Cu)
'
        'Line method:
'        m30_m = 0.269851
'        m30_b = 321.049
'        sum = OW + PL + Al + Cu
'
'        combined_SDE = m30_m * sum + m30_b
        
    Else 'All Others
        c0 = 1.48096
        c1 = -0.878005
        c2 = -0.716098
        c3 = 1.23318
        
'        combined_SDE = 1.25 * (c0 * OW + c1 * PL + c2 * Al + c3 * Cu)
        
        combined_SDE = (OW + PL + Al + Cu) * 0.0728819 + 404.345

        'Line method:
'        m30_m = 0.117222
'        m30_b = 362.638
'        sum = OW + PL + Al + Cu
'
'        combined_SDE = m30_m * sum + m30_b
    End If
    
    jboyd_calc_mixed_combined_SDE = combined_SDE
End Function

Public Function jboyd_calc_BmixPL_betaSDE(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c0, c1, c2, c3, c4, c5, c6, sde, Sum As Double
    
    Sum = OW + PL + Al + Cu
'
'    c0 = 428.18
'    c1 = -0.00521077
'    c2 = -0.0000213551
'    c3 = 0.00000000205956
'
'    sde = c0 + c1 * sum + c2 * sum * sum + c3 * sum * sum * sum
'
    c1 = 1.144
c2 = -0.977
c3 = -0.855
c4 = 5.395
    
    sde = c1 * OW + c2 * PL + c3 * Al + c4 * Cu
    
    jboyd_calc_BmixPL_betaSDE = sde
End Function

Function jboyd_mixed_beta_high_low_test(OW As Double, Al As Double)
    Dim m, b, lin_test As Double
    Dim beta_high As Integer
    
    m = 0.865203
    b = -183.48
    
    lin_test = m * OW + b
    
    If (Al > lin_test) Then
        beta_high = 0
    Else
        beta_high = 1
    End If
    
    jboyd_mixed_beta_high_low_test = beta_high
End Function

Public Function jboyd_calc_beta_rad_quality(OW As Double, PL As Double, Al As Double, Cu As Double, Optional return_type As Integer = 1, Optional input_pure_beta As Variant = -99, Optional rad_qual_mixed_input As Variant = -99)
    'Trying to differentiate between Sr and Kr-85 here.
    'Kr-85 is lower energy but has larger ratio values
    
    'Beta Low = 0
    'Beta High = 1
    
    'return_type = 1 --> Return an integer
    'return_type = 0 --> Return String
    
    Dim beta_low_high As Integer
    beta_low_high = 0
    Dim beta_low_high_string As String
    beta_low_high_string = ""
    
    Dim beta_indicator As Integer
    Dim pure_beta As Integer
    
    If (input_pure_beta = -99) Then
        pure_beta = jboyd_source_type_indicator(PL, Cu, "betas")
    Else 
        pure_beta = CInt(input_pure_beta)
    End If
    
    'Mixed sources
    Dim mixed_beta As Integer
    'If no input (mixed_beta_input) is provided, then we look up the mixed_beta condition.
    'mixed_beta = 0 --> NOT mixed
    'mixed_beta = 1 --> Mixed
    If (rad_qual_mixed_input = -99) Then
'        mixed_beta = jboyd_beta_indicator(OW, PL, Al, Cu, 1)
        mixed_beta = jboyd_source_type_indicator(PL, Cu, "mixed")
    Else
        mixed_beta = CInt(rad_qual_mixed_input)
    End If
    
    Dim signal_sum As Double
    
    Dim bh_sde, bl_sde As Double
    
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
    
'****************
'NOTE, WE CAN'T CLEANLY DISTINGUISH ALL LOW ENERGY BETA AND HIGH ENERGY BETA IN MIXED SOURCES
'There does seem to be a boundary between the two though in the OW element.
'We can draw a linear boundary between PL/Cu and OW and
'We can use a relation similar to: y = 0.057x - 28.5 (where y = PL/Cu and x = OW)
'If PL/Cu is larger than this value it will lie to the "left" of this line and most likely be a low-energy beta.
'Low energy betas will also exist above this but as of yet, there is no way to differentiate further

' BM = Beta Mixed (GENERAL)
' BLM = Low Energy Beta, Mixed
' BHM = High Energy Beta, Mixed

    Dim PL_Cu_linear_test As Double
    PL_Cu_linear_test = 0.057 * OW - 28.5
    
    If (OW_PL_SumScaled > 5) And (OW_Al_SumScaled > 5) And (OW_Cu_SumScaled > 5) Then
            beta_low_high_string = "BL"
    ElseIf ((PL_Cu >= 45) And (OW_Al_SumScaled > 0.2) And (OW_Cu > 45)) Or (OW_PL_CuCorr_SumScaled <= 0) Or (OW_PL > 20) Then
        beta_low_high = 0
        If (mixed_beta) And (pure_beta = 0) Then
            beta_low_high_string = "BLM"
'        ElseIf (OW_PL > 3) And (OW_Al > 3) And (OW_Cu > 3) Then

        Else
            beta_low_high_string = "BL"
        End If
    Else
        beta_low_high = 1
        If (mixed_beta) And (pure_beta = 0) Or ((PL_Cu >= 45) And (OW_Al_SumScaled > 0.2) And (OW_Cu > 35)) Or (OW_PL_CuCorr_SumScaled <= 0) Or (OW_PL > 20) And (pure_beta = 0) Then
            beta_low_high_string = "BHM"
        ElseIf (OW_Cu > 3) And (PL_Cu > 2) Then
            beta_low_high_string = "BH"
        Else
            beta_low_high_string = "Bmix"
        End If
    End If


    If (return_type) Then
        jboyd_calc_beta_rad_quality = beta_low_high
    Else
        jboyd_calc_beta_rad_quality = beta_low_high_string
    End If

End Function

Function jboyd_linear_beta_test(OW As Double, PL As Double, Cu As Double)
'   For the new PETG versions we can easily separate beta and no-betas using a linear relation between PL/Cu and OW
'   From the plot of PL/Cu vs OW, we can find a region separating data points from beta and non-beta sources.
'   The line is defined by y = 0.015x - 4.5 where y = PL/Cu and x = OW.
'   Everything below this line is "pure" and everything above is "beta"
'   This should be tested and stressed using much, much, more data.
    Dim PL_Cu, m, b, lin_test As Double
    PL_Cu = PL / Cu
    
    m = 0.015
    b = -4.5
    lin_test = m * OW + b
    
    If (PL_Cu > lin_test) Then
        jboyd_linear_beta_test = 1
    Else
        jboyd_linear_beta_test = 0
    End If
    
End Function

Public Function jboyd_calc_rad_quality(OW As Double, PL As Double, Al As Double, Cu As Double, Optional input_pure_beta As Variant = -99, _ 
        Optional input_mixed_beta As Variant = -99, Optional input_pure_photons As Variant = -99)
    
    Dim pure_beta, mixed_beta, pure_photons As Integer
    Dim rad_qual As String

    'New PL vs Cu tests
    If (input_pure_beta = -99) Then 
        pure_beta = jboyd_pure_beta_test(PL, Cu)
    Else 
        pure_beta = CInt(input_pure_beta)
    End If 

    If (input_mixed_beta = -99) Then 
        mixed_beta = jboyd_mixed_source_test(PL, Cu)
    Else 
        mixed_beta = CInt(input_mixed_beta)
    End If 

    If (input_pure_photons = -99) Then 
        pure_photons = jboyd_pure_photon_test(PL, Cu)
    Else 
        pure_photons = CInt(input_pure_photons)
    End If

    If (pure_beta = 1) Then
        rad_qual = jboyd_calc_beta_rad_quality(OW, PL, Al, Cu, 0, 0)
    End If
    
    If (pure_photons = 1) Then
        rad_qual = jboyd_calc_photon_rad_quality(OW, PL, Al, Cu, 0)
    End If
    
    If (mixed_beta = 1) Then
        rad_qual = "Bmix" + jboyd_calc_photon_rad_quality(OW, PL, Al, Cu, 0)
    End If
    
    jboyd_calc_rad_quality = rad_qual
        
End Function

Function jboyd_calc_mixed_DDE(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c0, c1, c2, c3, c4, c5, c6 As Double
    Dim mixed_dde As Double
        
c1 = -0.323
c2 = 0.497
c3 = 0.111
c4 = 0.604

    mixed_dde = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)
    
    jboyd_calc_mixed_DDE = mixed_dde
End Function

Function jboyd_calc_mixed_betaSDE(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c0, c1, c2, c3, c4, c5, c6 As Double
    Dim mixed_betaSDE As Double
        
c1 = 1.432
c2 = -1.676
c3 = 0.215
c4 = 0.226

    mixed_betaSDE = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)
    
    jboyd_calc_mixed_betaSDE = mixed_betaSDE
End Function

Public Function jboyd_calc_beta_DDE(OW As Double, PL As Double, Al As Double, Cu As Double, Optional input_beta_low_high As Variant = -99, Optional input_beta_mixed As Variant = -99)
    'beta_low_high = 99 -> Nothing specificed, will look it up
    'beta_low_high = 0 --> low energy beta
    'beta_low_high = 1 --> high energy beta
    Dim beta_low_high, beta_mixed As Integer

    If (beta_low_high = -99) Then
        beta_low_high = jboyd_calc_beta_rad_quality(OW, PL, Al, Cu, 1)
    Else 
        beta_low_high = CInt(input_beta_low_high)
    End If

    If (input_beta_mixed = -99) Then
        beta_mixed = jboyd_source_type_indicator(PL, Cu, "mixed")
    Else 
        beta_mixed = CInt(input_beta_mixed)
    End If
    
'    Can bake in a check that has the algorithm select the highest value of SDE for worst-case scenarious
    Dim BHM_DDE As Double
    Dim BH_DDE As Double
    Dim BLM_DDE As Double
    Dim BL_DDE As Double
    Dim BMgen_DDE As Double
    
    Dim beta_DDE As Double
    Dim c1 As Double
    Dim c2 As Double
    Dim c3 As Double
    Dim c4 As Double

    Dim pct As Double
    pct = 1#

    If (beta_low_high = 0) Then
        'low energy betas
        If (beta_mixed = 0) Then

            c1 = 0.000285893
            c2 = 0.00953786
            c3 = 0.00899577
            c4 = 0.979589
            BL_DDE = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)
        Else

            c1 = -0.217875
            c2 = 0.361691
            c3 = 0.0944204
            c4 = 0.633397
            
            BLM_DDE = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)
        End If
        
    ElseIf (beta_low_high = 1) Then
        'high energy betas
        If (beta_mixed = 0) Then

            c1 = -0.217875
            c2 = 0.361691
            c3 = 0.0944204
            c4 = 0.633397

            BH_DDE = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)
        Else
            c1 = -0.217875
            c2 = 0.361691
            c3 = 0.0944204
            c4 = 0.633397

            BHM_DDE = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)
        End If
        
    Else
        'General mixed beta case (77)
    c1 = -0.218
    c2 = 0.362
    c3 = 0.094
    c4 = 0.633
        BMgen_DDE = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)
    End If

    'Condition to check for pure beta --> Hard to identify pure from mixtures so we can look at total DDE calc and make some decisions
    'If the energies are 0 or less than 0 the we don't want that. We will just default it to ta standard low energy beta DDE val:
    
    jboyd_calc_beta_DDE = (Largest(beta_DDE, BLM_DDE, BHM_DDE, BL_DDE, BH_DDE, BMgen_DDE))
    
    If (jboyd_calc_beta_DDE <= 0) Then
    Dim Lc1, Lc2, Lc3, Lc4, Hc1, Hc2, Hc3, Hc4 As Double
    Dim gen_bl, gen_bh As Double
        Lc1 = 0.000285893
        Lc2 = 0.00953786
        Lc3 = 0.00899577
        Lc4 = 0.979589
        gen_bl = (Lc1 * OW) + (Lc2 * PL) + (Lc3 * Al) + (Lc4 * Cu)

        Hc1 = -0.217875
        Hc2 = 0.361691
        Hc3 = 0.0944204
        Hc4 = 0.633397
        gen_bh = (Hc1 * OW) + (Hc2 * PL) + (Hc3 * Al) + (Hc4 * Cu)
        
        jboyd_calc_beta_DDE = SmallestPositive(gen_bl, gen_h)
    End If

'    jboyd_calc_beta_DDE = (Largest(beta_DDE, BLM_DDE, BHM_DDE, BL_DDE, BH_DDE, BMgen_DDE)) * pct
    
End Function

Function jboyd_calc_BL_SDE(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c1, c2, c3, c4 As Double
    Dim Sum As Double
    Sum = OW + PL + Al + Cu
    
    'low energy betas
        'pure
       
'BASIC FUNCTION
'        c1 = 2.28162
'        c2 = -2.06136
'        c3 = 0#
'        c4 = 0.0116431
c1 = 1.535
c2 = 4.628
c3 = 25.375
c4 = -22.764

        bl_sde = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)
        
'ln+expo (Low energy beta - Kr ONLY)
'    c0 = 29.2469
'    c1 = 0.00005
'    c2 = 1.25
'    c3 = -1750
'    c4 = -600
'
'    bl_sde = c0 * Math.Log(sum / c1) + c2 * Exp((-0.5) * ((sum - c3) / c4))

'ln+expo (6 PAR) (Low energy beta - Kr ONLY)
'    c0 = 42.7603
'    c1 = 0.0000000476837
'    c2 = 0.005
'    c3 = 0.05
'    c4 = -3200
'    c5 = -487.5
'
'    If (sum <= 1.001 * c1) Then
'        sum = 1.00105 * sum
'    End If
'
'    bl_sde = c0 * Math.Log((sum - c1) / c2) + c3 * Exp((-0.5) * ((sum - c4) / c5))
        jboyd_calc_BL_SDE = bl_sde
End Function

Function jboyd_calc_BLM_SDE(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c1, c2, c3, c4 As Double
    Dim Sum As Double
    Sum = OW + PL + Al + Cu
    
    'low energy betas
        'mixed
        c0 = 611
        c1 = 0.000000123998
        c2 = 0.358
        c3 = 1304
        
        BLM_SDE = c0 + (c1) * ((c2 * Sum - c3) ^ (3))
        
        jboyd_calc_BLM_SDE = BLM_SDE
End Function

Function jboyd_calc_BH_SDE(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c1, c2, c3, c4 As Double
    Dim Sum As Double
    Sum = OW + PL + Al + Cu

    'high energy betas
        'pure

'BASIC FUNCTION
'        c1 = 1.5
'        c2 = -1.5
'        c3 = 0.3
'        c4 = -0.3
'
'        bh_sde = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)

''Irreducial Form Factor Pol3 (high energy beta - Sy ONLY)
'    c0 = -649.878
'    c1 = 0.0000000180641
'    c2 = -0.000135992
'    c3 = 0.40616
'    bh_sde = (sum - c0) * (c1 * sum * sum + c2 * sum + c3)

'ln+expo (high energy beta - Sy ONLY)
    c0 = 27.9222
    c1 = 0.0001
    c2 = 1.5
    c3 = -2500
    c4 = -664.077

    bh_sde = c0 * Math.Log(Sum / c1) + c2 * Exp((-0.5) * ((Sum - c3) / c4))

    jboyd_calc_BH_SDE = bh_sde

End Function

Function jboyd_calc_BHM_SDE(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c1, c2, c3, c4 As Double
    Dim Sum As Double
    Sum = OW + PL + Al + Cu
    
    'high energy betas
        'mixed
'        c1 = 0.76
'        c2 = -0.03
'        c3 = -0.55
'        c4 = 1.06
'        BHM_SDE = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)
        
'ln+expo (high energy beta - Sy ONLY)
    c0 = 29.2469
    c1 = 0.0001
    c2 = 0.00005
    c3 = -1750
    c4 = -600

    BHM_SDE = c0 * Math.Log(Sum / c1) + c2 * Exp((-0.5) * ((Sum - c3) / c4))
    
        jboyd_calc_BHM_SDE = BHM_SDE
End Function

Function jboyd_calc_BM_SDE(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c1, c2, c3, c4 As Double
    Dim Sum As Double
    Sum = OW + PL + Al + Cu
End Function

Function jboyd_pure_beta_test(PL As Double, Cu As Double)
    Dim m1, m2, b1, b2, upper_test, lower_test As Double
    Dim pure_beta As Integer
    
    'Based on a plot of PL vs Cu for photons, betas, and mixed, we can find some boundaries between each of these types
    'We can isolate most of the pure betas between two lines. Values between the lines are pure beta.
    'upper:
    m1 = 9.2329
    b1 = 17.5
    upper_test = m1 * Cu + b1
    
    'lower:
    m2 = 9.2329
    b2 = -250
    lower_test = m2 * 9.2329 + b2
    
    'I will also put a cap on Cu just so that we "truncate" the region between the two lines
    'If it goes to infinity it may overlap with the region defined by the pol4 for pure photons
    
    If (PL < upper_test) And (PL > lower_test) And (Cu < 200) Then
        pure_beta = 1
    Else
        pure_beta = 0
    End If
    
    jboyd_pure_beta_test = pure_beta
    
End Function

Function jboyd_pure_photon_test(PL As Double, Cu As Double)
    Dim c0, c1, c2, c3, pol4_test As Double
    Dim pure_photon As Integer

    'Based on a plot of PL vs Cu for photons, betas, and mixed, we can find some boundaries between each of these types
    'The boundary between pure and (mixed or beta) can be defined by a 4th order polynomial.
    '->Values below the pol4 are mixed or beta, and values above the pol4 are pure photons
    ' y = A + B(x-D) + C(x-D)^2, where y = PL and x = Cu
    'Pol4:
    c0 = 311.6
    c1 = -2.6
    c2 = 0.0183
    c3 = 534
    
    pol4_test = c0 + c1 * (Cu - c3) + c2 * ((Cu - c3) ^ (2))
    
    If (PL > pol4_test) Then
        pure_photon = 1
    Else
        pure_photon = 0
    End If
    
    jboyd_pure_photon_test = pure_photon
End Function

Function jboyd_mixed_source_test(PL As Double, Cu As Double)
    Dim pure_photon, pure_beta, mixed_source As Integer
    
    pure_photon = jboyd_pure_photon_test(PL, Cu)
    pure_beta = jboyd_pure_beta_test(PL, Cu)
    
    If (pure_photon = 0) And (pure_beta = 0) Then
        mixed_source = 1
    Else
        mixed_source = 0
    End If
    
    jboyd_mixed_source_test = mixed_source
End Function

Function jboyd_source_type_indicator(PL As Double, Cu As Double, Optional return_type As String = "general")
    Dim beta_indicator, pure_beta, mixed_beta, pure_photon As Integer
    Dim pure_beta_string, mixed_beta_string, pure_photon_string As String
    Dim source_type As String
    
    pure_beta = jboyd_pure_beta_test(PL, Cu)
    pure_photon = jboyd_pure_photon_test(PL, Cu)
    mixed_beta = jboyd_mixed_source_test(PL, Cu)
    
    If (pure_beta = 1) Or (mixed_beta = 1) Then
        beta_indicator = 1
    Else
        beta_indicator = 0
    End If
    
    If (pure_beta = 1) Then
        source_type = "pure_betas"
    End If
    If (pure_photon = 1) Then
        source_type = "pure_photons"
    End If
    If (mixed_beta = 1) Then
        source_type = "mixed_source"
    End If
    
    If (StrComp("beta_indicator", return_type, vbBinaryCompare) = 0) Then
        jboyd_source_type_indicator = beta_indicator
    ElseIf (StrComp("photons", return_type, vbBinaryCompare) = 0) Then
        jboyd_source_type_indicator = pure_photon
    ElseIf (StrComp("betas", return_type, vbBinaryCompare) = 0) Then
        jboyd_source_type_indicator = pure_beta
    ElseIf (StrComp("mixed", return_type, vbBinaryCompare) = 0) Then
        jboyd_source_type_indicator = mixed_beta
    Else
        jboyd_source_type_indicator = source_type
    End If
    

End Function

Function jboyd_calc_mixed_beta_SDE(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c1, c2, c3, c4, mixed_beta_SDE As Double
    
'Calibrated to PETG 0.04 data
c1 = 1.529
c2 = -1.736
c3 = 0.118
c4 = 0.322

        mixed_beta_SDE = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)

    
        jboyd_calc_mixed_beta_SDE = mixed_beta_SDE
End Function

Function jboyd_calc_mixed_beta_SDE_PL(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c0, c1, c2, c3, c4, c5, c6, c7, mixed_beta_SDE, Sum As Double
    Sum = OW + PL + Al + Cu
    
''Calibrated to PETG 0.04 data
c0 = 304.84
c1 = 0.352801
c2 = -0.000409837
c3 = 0.000000196393
c4 = -5.05617E-11
c5 = 7.01519E-15
c6 = -4.84321E-19
c7 = 1.27161E-23

        mixed_beta_SDE = c0 + c1 * (Sum) + c2 * ((Sum) ^ (2)) + c3 * ((Sum) ^ (3)) + c4 * ((Sum) ^ (4)) + c5 * ((Sum) ^ (5)) + c6 * ((Sum) ^ (6)) + c7 * ((Sum) ^ (7))
'
'c1 = 1.179
'c2 = -0.973
'c3 = -0.936
'c4 = 5.498
'    mixed_beta_SDE = c1 * OW + c2 * PL + c3 * Al + c4 * Cu

        jboyd_calc_mixed_beta_SDE_PL = mixed_beta_SDE
End Function

Function jboyd_calc_mixed_beta_SDE_PH(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c0, c1, c2, c3, c4, c5, c6, c7, mixed_beta_SDE, Sum As Double
    Sum = OW + PL + Al + Cu
    
''Calibrated to PETG 0.04 data
'c0 = 513.812
'c1 = -0.244742
'c2 = -2.26115e-007
'
'    mixed_beta_SDE = c0 + c1 * sum + c2 * (sum * sum)

c1 = 2.28
c2 = -1.273
c3 = 0.004
c4 = -1.022
    mixed_beta_SDE = c1 * OW + c2 * PL + c3 * Al + c4 * Cu

        jboyd_calc_mixed_beta_SDE_PH = mixed_beta_SDE
End Function

Function jboyd_calc_BL_beta_combined_SDE(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c1, c2, c3, c4, Sum, BL_beta_SDE As Double
    Dim m, b As Double
    
    'Can try two methods for calculating the combined SDE for BL betas:
    'Standard Coefficient method:
    'SEEMS TO BE THE MORE ACCURATE METHOD
c1 = 2.366
c2 = -3.125
c3 = -0.005
c4 = 0.941

    BL_beta_SDE = c1 * OW + c2 * PL + c3 * Al + c4 * Cu
    
    'Linear fit to "Combined SDE vs (OW + PL + Al + Cu)":
    'Sum:
'    Sum = OW + PL + Al + Cu
'    'Slope:
'    m = 2.407
'    'Intercept:
'    b = -159.022
'
'    BL_beta_SDE = m * Sum + b

    jboyd_calc_BL_beta_combined_SDE = BL_beta_SDE
End Function

Function jboyd_calc_BH_beta_combined_SDE(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c1, c2, c3, c4, BH_beta_SDE As Double

    
    'Standard Coefficient method:

    c1 = 1.487
    c2 = -1.534
    c3 = 0.318
    c4 = 0.658
    
    BH_beta_SDE = c1 * OW + c2 * PL + c3 * Al + c4 * Cu
    
    jboyd_calc_BH_beta_combined_SDE = BH_beta_SDE
End Function

Function jboyd_calc_combined_SDE(OW As Double, PL As Double, Al As Double, Cu As Double, _
    Optional input_beta_indicator As Variant = -99, Optional input_pure_beta As Variant = -99, Optional input_beta_mixed As Variant = -99, _
    Optional input_pure_photons As Variant = -99, Optional input_beta_low_high As Variant = -99)
    
    Dim combined_SDE As Double
    Dim beta_indicator, pure_beta, beta_mixed, pure_photons, beta_low_high As Integer
    
    'If user-inputs for the indicators are NOT provided then we will calculate them so that we can pass them along
    If (input_beta_indicator = -99) Then
        beta_indicator = jboyd_source_type_indicator(PL, Cu, "beta_indicator")
    Else 
        beta_indicator = CInt(input_beta_indicator)
    End If
    
    If (input_pure_beta = -99) Then
        pure_beta = jboyd_source_type_indicator(PL, Cu, "betas")
    Else 
        pure_beta = CInt(input_pure_beta)
    End If

    If (input_beta_mixed = -99) Then
        beta_mixed = jboyd_source_type_indicator(PL, Cu, "mixed")
    Else 
        beta_mixed = CInt(input_beta_mixed)
    End If

    If (input_pure_photons = -99) Then
        pure_photons = jboyd_source_type_indicator(PL, Cu, "photons")
    Else 
        pure_photons = CInt(input_pure_photons)
    End If

    If (input_beta_low_high = -99) Then
        beta_low_high = jboyd_calc_beta_rad_quality(OW, PL, Al, Cu, 1)
    Else 
        beta_low_high = CInt(input_beta_low_high)
    End If
    
    'Pure Photons are simplest case (No sub-divisions of photons for the SDE calculations)
    If (pure_photons = 1) Then
        combined_SDE = jboyd_calc_photons_combined_SDE(OW, PL, Al, Cu)
    Else
        'So we know have betas/mixed
        'Mixed sources don't have any sub-divisions for the combined SDE calculations:
        'I am not doing If/Else statements here because I want to explicitly declare my conditions...
        'just a preferential decision here....

        If (beta_mixed = 1) Then
            combined_SDE = jboyd_calc_mixed_combined_SDE(OW, PL, Al, Cu)
        End If
        If (pure_beta = 1) Then
            'We have pure betas but we need to classify the calculation based on high (BH) or low (BL) energy:
            If (beta_low_high = 0) Then 'BL:
                combined_SDE = jboyd_calc_BL_beta_combined_SDE(OW, PL, Al, Cu)
            End If
            If (beta_low_high = 1) Then 'BH:
                combined_SDE = jboyd_calc_BH_beta_combined_SDE(OW, PL, Al, Cu)
            End If
        End If
    End If
    
    jboyd_calc_combined_SDE = combined_SDE
    
End Function

Function jboyd_calc_beta_SDE(OW As Double, PL As Double, Al As Double, Cu As Double, Optional beta_low_high As Integer = -99, Optional input_beta_indicator As Variant = -99, _ 
            Optional input_pure_beta As Variant = -99, Optional input_pure_photons As Variant = -99)
    'beta_low_high = 99 -> Nothing specificed, will look it up
    'beta_low_high = 0 --> low energy beta
    'beta_low_high = 1 --> high energy beta
    'beta_low_high = 77 --> Beta mixed (General)
    
    If (beta_low_high < 0) Then
        beta_low_high = jboyd_calc_beta_rad_quality(OW, PL, Al, Cu, 1)
    End If
    
    If (input_beta_indicator = -99) Then
        beta_indicator = jboyd_source_type_indicator(PL, Cu, "beta_indicator")
    Else 
        beta_indicator = CInt(input_beta_indicator)
    End If 

    If (input_pure_beta = -99) Then
        pure_beta = jboyd_source_type_indicator(PL, Cu, "betas")
    Else
        pure_beta = CInt(input_pure_beta)
    End If 

    If (input_pure_photons = -99) Then
        pure_photons = jboyd_source_type_indicator(PL, Cu, "photons")
    Else
        pure_photons = CInt(input_pure_photons)
    End If 

    If (beta_low_high >= 0) Then
        beta_mixed = 1
    Else
        beta_mixed = jboyd_mixed_source_test(PL, Al)
    End If
    
    'Can bake in a check that has the algorithm select the highest value of SDE for worst-case scenarious
    Dim BHM_SDE As Double
    Dim bh_sde As Double
    Dim BLM_SDE As Double
    Dim bl_sde As Double
    
    Dim final_MAX_beta_SDE As Double

    bh_sde = jboyd_calc_BH_SDE(OW, PL, Al, Cu)
    bl_sde = jboyd_calc_BL_SDE(OW, PL, Al, Cu)
    BHM_SDE = jboyd_calc_BHM_SDE(OW, PL, Al, Cu)
    BLM_SDE = jboyd_calc_BLM_SDE(OW, PL, Al, Cu)

    If (pure_photons = 1) Then
        jboyd_calc_beta_SDE = 0
    Else
        If (beta_low_high = 0) Then
            'low energy betas
    '        If (beta_mixed = 0) Then
            If (pure_beta = 1) Then
                final_MAX_beta_SDE = bl_sde
            Else
                final_MAX_beta_SDE = BLM_SDE
            End If
        Else
            'high energy betas
    '        If (beta_mixed = 0) Then
            If (pure_beta = 1) Then
                final_MAX_beta_SDE = bh_sde
            Else
                final_MAX_beta_SDE = BHM_SDE
            End If
        End If
    End If
    
'    jboyd_calc_beta_SDE = Largest(bl_sde, BLM_SDE, bh_sde, BHM_SDE)
    jboyd_calc_beta_SDE = final_MAX_beta_SDE
    
End Function

Function jboyd_calc_photon_SDE_for_BL(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c0, c1, c2, c3, c4, c5, c6 As Double
        
c1 = 0.003
c2 = 0.398
c3 = 0.088
c4 = 0.492

    ph_bl_sde = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)
    
    jboyd_calc_photon_SDE_for_BL = ph_bl_sde
End Function

Function jboyd_calc_photon_SDE_for_BH(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c0, c1, c2, c3, c4, c5, c6 As Double
        
c1 = 0.023
c2 = -0.028
c3 = -0.009
c4 = 0.936

    ph_bl_sde = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)
    
    jboyd_calc_photon_SDE_for_BH = ph_bl_sde
End Function

Function jboyd_calc_photon_SDE_for_purePH(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c0, c1, c2, c3, c4, c5, c6 As Double
        
c1 = 0.406
c2 = 0.764
c3 = -1.47
c4 = 1.285

    ph_bl_sde = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)
    
    jboyd_calc_photon_SDE_for_purePH = ph_bl_sde
End Function

Function jboyd_calc_photon_SDE_for_lowE_purePH(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c0, c1, c2, c3, c4, c5, c6 As Double

c1 = 0.14
c2 = 0.696
c3 = -0.379
c4 = -4.743

    ph_bl_sde = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)
    
    jboyd_calc_photon_SDE_for_lowE_purePH = ph_bl_sde
End Function

Function jboyd_calc_photon_SDE_for_highE_purePH(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c0, c1, c2, c3, c4, c5, c6 As Double

c1 = 0.023
c2 = -0.028
c3 = -0.009
c4 = 0.936

    ph_bl_sde = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)
    
    jboyd_calc_photon_SDE_for_lowE_purePH = ph_bl_sde
End Function

Public Function jboyd_calc_photon_Cu_SDE(PL As Double, Al As Double, Cu As Double)
    Dim PL_Cu As Double
    Dim Al_Cu As Double
    Dim photon_Cu_SDE As Double
    Dim Fact As Double
    
    Dim c0, c1, c2, c3, c4, c5, c6, pct As Double
    
c0 = -0.067041
c1 = 0.0035531
c2 = -0.000041356
c3 = 0.852763
c4 = 0.0404893
c5 = 0.0149254
    
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
    
    jboyd_calc_photon_Cu_SDE = photon_Cu_SDE
End Function

Public Function jboyd_calc_photon_Cu_SDE_purePH(PL As Double, Al As Double, Cu As Double)
    Dim PL_Cu As Double
    Dim Al_Cu As Double
    Dim photon_Cu_SDE As Double
    Dim Fact As Double
    
    Dim c0, c1, c2, c3, c4, c5, c6, pct As Double
    
c0 = -0.067041
c1 = 0.0035531
c2 = -0.000041356
c3 = 0.852763
c4 = 0.0404893
c5 = 0.0149254
    
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
    
    jboyd_calc_photon_Cu_SDE_purePH = photon_Cu_SDE
End Function

Public Function calc_photon_PL_SDE(PL As Double, Cu As Double)

    Dim PL_Cu As Double
    Dim photon_PL_SDE As Double

    PL_Cu = PL / Cu

    Dim c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10 As Double

    c0 = -0.388614
    c1 = 6.63087
    c2 = -0.15941
    c3 = 0.407571
    c4 = -0.154127
    c5 = 0.0404295
    c6 = 7.17117
    c7 = 0.892884
    c8 = -2.48003
    c9 = 0.599806
    c10 = 0.0413942
    
    photon_PL_SDE = ((c0 + _
                (c1 * Math.Log(PL_Cu)) + _
                (c2 * (Math.Log(PL_Cu) ^ 2)) + _
                (c3 * (Math.Log(PL_Cu) ^ 3)) + _
                (c4 * (Math.Log(PL_Cu) ^ 4)) + _
                (c5 * (Math.Log(PL_Cu) ^ 5))) / _
                (1# + (c6 * Math.Log(PL_Cu)) + _
                (c7 * (Math.Log(PL_Cu) ^ 2)) + _
                (c8 * (Math.Log(PL_Cu) ^ 3)) + _
                (c9 * (Math.Log(PL_Cu) ^ 4)) + _
                (c10 * (Math.Log(PL_Cu) ^ 5)))) * PL
                
    calc_photon_PL_SDE = photon_PL_SDE
End Function

Public Function jboyd_calc_photon_PL_SDE(PL As Double, Cu As Double)

    Dim PL_Cu As Double
    Dim photon_PL_SDE As Double

    PL_Cu = PL / Cu
    
    Dim c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, pct As Double

c0 = 1.03317
c1 = 0.899083
c2 = 0.0414097
c3 = 0.128686
c4 = -0.165167
c5 = 0.0285842
c6 = 1.25184
c7 = 5.04014
c8 = -2.929
c9 = 0.396408
c10 = 0.000304115

    photon_PL_SDE = ((c0 + _
                (c1 * Math.Log(PL_Cu)) + _
                (c2 * (Math.Log(PL_Cu) ^ 2)) + _
                (c3 * (Math.Log(PL_Cu) ^ 3)) + _
                (c4 * (Math.Log(PL_Cu) ^ 4)) + _
                (c5 * (Math.Log(PL_Cu) ^ 5))) / _
                (1# + (c6 * Math.Log(PL_Cu)) + _
                (c7 * (Math.Log(PL_Cu) ^ 2)) + _
                (c8 * (Math.Log(PL_Cu) ^ 3)) + _
                (c9 * (Math.Log(PL_Cu) ^ 4)) + _
                (c10 * (Math.Log(PL_Cu) ^ 5)))) * PL
                
                
    jboyd_calc_photon_PL_SDE = photon_PL_SDE
End Function

Public Function jboyd_calc_photon_PL_SDE_purePH(PL As Double, Cu As Double)

    Dim PL_Cu As Double
    Dim photon_PL_SDE As Double

    PL_Cu = PL / Cu
    
    Dim c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, pct As Double

c0 = 1.03317
c1 = 0.899083
c2 = 0.0414097
c3 = 0.128686
c4 = -0.165167
c5 = 0.0285842
c6 = 1.25184
c7 = 5.04014
c8 = -2.929
c9 = 0.396408
c10 = 0.000304115

    photon_PL_SDE = ((c0 + _
                (c1 * Math.Log(PL_Cu)) + _
                (c2 * (Math.Log(PL_Cu) ^ 2)) + _
                (c3 * (Math.Log(PL_Cu) ^ 3)) + _
                (c4 * (Math.Log(PL_Cu) ^ 4)) + _
                (c5 * (Math.Log(PL_Cu) ^ 5))) / _
                (1# + (c6 * Math.Log(PL_Cu)) + _
                (c7 * (Math.Log(PL_Cu) ^ 2)) + _
                (c8 * (Math.Log(PL_Cu) ^ 3)) + _
                (c9 * (Math.Log(PL_Cu) ^ 4)) + _
                (c10 * (Math.Log(PL_Cu) ^ 5)))) * PL
                
                
    jboyd_calc_photon_PL_SDE_purePH = photon_PL_SDE
End Function

Public Function jboyd_calc_average_photon_SDE(OW As Double, PL As Double, Al As Double, Cu As Double, Optional input_gen_rad_qual As Variant = -99, _ 
        Optional input_photon_rad_qual As Variant = -99, Optional input_pure_beta As Variant = -99, Optional input_mixed_beta As Variant = -99, _ 
            Optional input_pure_photon As Variant = -99, Optional input_high_energy_beta As Variant = -99)
    Dim c1, c2, c3, c4 As Double
    Dim photon_Cu, photon_Cu_purePH As Double
    Dim photon_PL, photon_PL_purePH As Double
    Dim average_photon_SDE, average_photon_SDE_purePH As Double
    Dim pure_beta, mixed_beta, pure_photon, general_beta_indicator As Integer
    Dim high_energy_beta As Integer
    Dim general_rad_qual, photon_rad_qual As String
    
    If (input_photon_rad_qual = -99) Then
        photon_rad_qual = jboyd_calc_photon_rad_quality(OW, PL, Al, Cu, 0)
    Else 
        photon_rad_qual = input_photon_rad_qual
    End If 

    If (input_pure_beta = -99) Then
        pure_beta = jboyd_pure_beta_test(PL, Cu)
    Else 
        pure_beta = CInt(input_pure_beta)
    End If 

    If (input_mixed_beta = -99) Then
        mixed_beta = jboyd_mixed_source_test(PL, Cu)
    Else 
        mixed_beta = CInt(input_mixed_beta)
    End If 

    If (input_pure_photon = -99) Then
        pure_photon = jboyd_pure_photon_test(PL, Cu)
    Else 
        pure_photon = CInt(pure_photon)
    End If 

    If (input_high_energy_beta = -99) Then
        high_energy_beta = jboyd_mixed_beta_high_low_test(OW, Al) '0 if BL and 1 if BH
    Else
        high_energy_beta = CInt(input_high_energy_beta)
    End If
    
    If (pure_beta = 1) Or (mixed_beta = 1) Then
        general_beta_indicator = 1
    Else
        general_beta_indicator = 0
    End If

    If (input_gen_rad_qual = -99) Then
        general_rad_qual = jboyd_calc_rad_quality(OW, PL, Al, Cu, pure_beta, mixed_beta, pure_photon)
    Else 
        general_rad_qual = input_gen_rad_qual
    End If
    
    photon_Cu = jboyd_calc_photon_Cu_SDE(PL, Al, Cu)
    photon_PL = jboyd_calc_photon_PL_SDE(PL, Cu)
    
    photon_Cu_purePH = jboyd_calc_photon_Cu_SDE_purePH(PL, Al, Cu)
    photon_PL_purePH = jboyd_calc_photon_PL_SDE_purePH(PL, Cu)
    
    average_photon_SDE = (photon_Cu + photon_PL) / 2#
    average_photon_SDE_purePH = (photon_Cu_purePH + photon_PL_purePH) / 2
    
    
    If (general_beta_indicator = 1) Then 'Betas exist: either pure or mixed
        If (pure_beta = 1) Then
            If (high_energy_beta = 0) Then  'Low energy beta
                jboyd_calc_average_photon_SDE = 1.05 * jboyd_calc_photon_SDE_for_BL(OW, PL, Al, Cu)
            Else 'high energy beta
                jboyd_calc_average_photon_SDE = jboyd_calc_photon_SDE_for_BH(OW, PL, Al, Cu)
            End If
        Else 'mixed beta
'            jboyd_calc_average_photon_SDE = jboyd_calc_mixed_beta_SDE(OW, PL, Al, Cu)
            c1 = -0.045
            c2 = 1.312
            c3 = -1.568
            c4 = 1.255

            jboyd_calc_average_photon_SDE = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)
            
        End If
        
    Else 'pure photons only
'        If (StrComp("PL", photon_rad_qual, vbBinaryCompare) = 0) And (mixed_beta <> 1) Then
'            jboyd_calc_average_photon_SDE = jboyd_calc_photon_SDE_for_lowE_purePH(OW, PL, Al, Cu)
'        ElseIf (mixed_beta = 1) Then
'            jboyd_calc_average_photon_SDE = jboyd_calc_mixed_beta_SDE(OW, PL, Al, Cu)
'        Else
'            jboyd_calc_average_photon_SDE = jboyd_calc_photon_SDE_for_purePH(OW, PL, Al, Cu)
'        End If
            c1 = -0.045
            c2 = 1.312
            c3 = -1.568
            c4 = 1.255

            jboyd_calc_average_photon_SDE = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)
    End If


End Function

Function calc_mixed_DDE_gaus(OW As Double, PL As Double, Cu As Double)
    Dim OWCu_PLCu_ratio As Double
    Dim A, mu, sigma, offset, gaus As Double
    
    If (PL = Cu) Then
        PL = 1.0000001 * PL
    End If
    
    OWCu_PLCu_ratio = (OW - Cu) / (PL - Cu)
    
    A = 1857
    mu = 1
    sigma = 0.102
    offset = 512
    
    gaus = A * Exp((-0.5) * (((OWCu_PLCu_ratio - mu) / sigma) ^ (2))) + offset
    calc_mixed_DDE_gaus = gaus
End Function


Public Function calc_average_photon_SDE(PL As Double, Al As Double, Cu As Double)
    Dim photon_Cu As Double
    Dim photon_PL As Double
    Dim average_photon_SDE As Double
    
    photon_Cu = calc_photon_Cu_SDE(PL, Al, Cu)
    photon_PL = calc_photon_PL_SDE(PL, Cu)
    
    average_photon_SDE = (photon_Cu + photon_PL) / 2#
    
    calc_average_photon_SDE = average_photon_SDE
    
End Function

Public Function calc_photon_Cu_DDE(Al As Double, Cu As Double)
    Dim Al_Cu, c0, c1 As Double
    Dim photon_Cu_DDE As Double
    
    Al_Cu = Al / Cu
    
    c0 = 0.729750759967566
    c1 = 0.129198337923333
    
    photon_Cu_DDE = ((Abs((c0 + (c1 * ((dblAlCu) ^ 2))))) ^ (0.5)) * Cu
    
    calc_photon_Cu_DDE = photon_Cu_DDE

End Function

Public Function jboyd_calc_photon_Cu_DDE_energyPL(Al As Double, Cu As Double)
    Dim Al_Cu, c0, c1 As Double
    Dim photon_Cu_DDE_energyPL As Double
    
    Al_Cu = Al / Cu
'
'c0 = 1.46086
'c1 = 4.51799
    
c0 = 1.31508
c1 = 4.08069
    photon_Cu_DDE_energyPL = ((Abs((c0 + (c1 * ((dblAlCu) ^ 2))))) ^ (0.5)) * Cu
    
    jboyd_calc_photon_Cu_DDE_energyPL = 1.15 * photon_Cu_DDE_energyPL

End Function

Public Function jboyd_calc_photon_Cu_DDE_energyPM(Al As Double, Cu As Double)
    Dim Al_Cu, c0, c1 As Double
    Dim photon_Cu_DDE_energyPM As Double
    
    Al_Cu = Al / Cu
    
c0 = 0.000192437
c1 = 0.000026903
    
    photon_Cu_DDE_energyPM = ((Abs((c0 + (c1 * ((dblAlCu) ^ 2))))) ^ (0.5)) * Cu
    
    jboyd_calc_photon_Cu_DDE_energyPM = photon_Cu_DDE_energyPM

End Function

Public Function jboyd_calc_photon_Cu_DDE_energyPH(Al As Double, Cu As Double)
    Dim Al_Cu, c0, c1 As Double
    Dim photon_Cu_DDE_energyPH As Double
    
    Al_Cu = Al / Cu
    
    c0 = 0.0000144904
    c1 = 1.99832
    
    photon_Cu_DDE_energyPH = ((Abs((c0 + (c1 * ((dblAlCu) ^ 2))))) ^ (0.5)) * Cu
    
    jboyd_calc_photon_Cu_DDE_energyPH = photon_Cu_DDE_energyPH

End Function

Public Function calc_photon_PL_DDE(PL As Double, Cu As Double)

    Dim PL_Cu As Double
    Dim photon_PL_DDE As Double

    PL_Cu = PL / Cu
    
    Dim c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, pct As Double
    pct = 1#
    c0 = 2.27875800359032 * pct
    c1 = -0.160519484455387 * pct
    c2 = 0.444931311949578 * pct
    c3 = -1.39776305022687E-02 * pct
    c4 = 1.25778369690611E-04 * pct
    c5 = -6.65137516928388E-02 * pct
    c6 = 1.49910566096546 * pct
    c7 = -4.14191308266288E-02 * pct
    c8 = 2.2610007373255E-04 * pct
    c9 = 2.1999287177623E-06 * pct
    
    photon_PL_DDE = ((c0 + _
                (c1 * (PL_Cu)) + _
                (c2 * ((PL_Cu) ^ 2)) + _
                (c3 * ((PL_Cu) ^ 3)) + _
                (c4 * ((PL_Cu) ^ 4))) / _
                (1# + (c5 * (PL_Cu)) + _
                (c6 * ((PL_Cu) ^ 2)) + _
                (c7 * ((PL_Cu) ^ 3)) + _
                (c8 * ((PL_Cu) ^ 4)) + _
                (c9 * ((PL_Cu) ^ 5)))) * PL
                
                
    calc_photon_PL_DDE = photon_PL_DDE
End Function

Public Function jboyd_calc_photon_PL_DDE_energyPL(PL As Double, Cu As Double)

    Dim PL_Cu As Double
    Dim photon_PL_DDE_energyPL As Double

    PL_Cu = PL / Cu
    
    Dim c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, pct As Double

    c0 = 1.09282
    c1 = -1.59537
    c2 = 0.462245
    c3 = -0.00386405
    c4 = 0.000424096
    c5 = -1.61565
    c6 = 0.580476
    c7 = 0.0232276
    c8 = 0.000426323
    c9 = 0.00001267
    
    photon_PL_DDE_energyPL = ((c0 + _
                (c1 * (PL_Cu)) + _
                (c2 * ((PL_Cu) ^ 2)) + _
                (c3 * ((PL_Cu) ^ 3)) + _
                (c4 * ((PL_Cu) ^ 4))) / _
                (1# + (c5 * (PL_Cu)) + _
                (c6 * ((PL_Cu) ^ 2)) + _
                (c7 * ((PL_Cu) ^ 3)) + _
                (c8 * ((PL_Cu) ^ 4)) + _
                (c9 * ((PL_Cu) ^ 5)))) * PL
                
                
    jboyd_calc_photon_PL_DDE_energyPL = 1.15 * photon_PL_DDE_energyPL
End Function

Public Function jboyd_calc_photon_PL_DDE_energyPM(PL As Double, Cu As Double)

    Dim PL_Cu As Double
    Dim photon_PL_DDE_energyPM As Double

    PL_Cu = PL / Cu
    
    Dim c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, pct As Double

    c0 = 1.59703
    c1 = 0.437369
    c2 = 0.872299
    c3 = -0.0170836
    c4 = 0.000136882
    c5 = -1.21716
    c6 = 1.73392
    c7 = -0.0323779
    c8 = 0.000361399
    c9 = -0.000000199984
    
    photon_PL_DDE_energyPM = ((c0 + _
                (c1 * (PL_Cu)) + _
                (c2 * ((PL_Cu) ^ 2)) + _
                (c3 * ((PL_Cu) ^ 3)) + _
                (c4 * ((PL_Cu) ^ 4))) / _
                (1# + (c5 * (PL_Cu)) + _
                (c6 * ((PL_Cu) ^ 2)) + _
                (c7 * ((PL_Cu) ^ 3)) + _
                (c8 * ((PL_Cu) ^ 4)) + _
                (c9 * ((PL_Cu) ^ 5)))) * PL
                
                
    jboyd_calc_photon_PL_DDE_energyPM = photon_PL_DDE_energyPM
End Function

Public Function jboyd_calc_photon_PL_DDE_energyPH(PL As Double, Cu As Double)

    Dim PL_Cu As Double
    Dim photon_PL_DDE_energyPH As Double

    PL_Cu = PL / Cu
    
    Dim c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, pct As Double
    c0 = 0.322244
    c1 = -2.33279
    c2 = -4.5112
    c3 = 4.2319
    c4 = 2.84632
    c5 = -3.01483
    c6 = 1.31702
    c7 = -0.584928
    c8 = -0.922422
    c9 = 2.48616
    
    photon_PL_DDE_energyPH = ((c0 + _
                (c1 * (PL_Cu)) + _
                (c2 * ((PL_Cu) ^ 2)) + _
                (c3 * ((PL_Cu) ^ 3)) + _
                (c4 * ((PL_Cu) ^ 4))) / _
                (1# + (c5 * (PL_Cu)) + _
                (c6 * ((PL_Cu) ^ 2)) + _
                (c7 * ((PL_Cu) ^ 3)) + _
                (c8 * ((PL_Cu) ^ 4)) + _
                (c9 * ((PL_Cu) ^ 5)))) * PL
                
                
    jboyd_calc_photon_PL_DDE_energyPH = photon_PL_DDE_energyPH
End Function

Public Function calc_average_photon_DDE(PL As Double, Al As Double, Cu As Double, Optional photon_energy As String = "all")
    'Options for photon_energy are:
    ' all : coefficients calibrated from datasets considering ALL energies, i.e., PL, PM, PH
    ' PL : coefficients calibrated from datasets only with PL
    ' PM ""
    ' PH ""
    Dim photon_Cu As Double
    Dim photon_PL As Double
    Dim average_photon_DDE As Double
    
    photon_Cu = calc_photon_Cu_DDE(Al, Cu)
    photon_PL = calc_photon_PL_DDE(PL, Cu)

    average_photon_DDE = (photon_Cu + photon_PL) / 2#
    
    calc_average_photon_DDE = average_photon_DDE
    
End Function

Public Function jboyd_calc_photon_Cu_DDE(Al As Double, Cu As Double)
    Dim Al_Cu, c0, c1 As Double
    Dim photon_Cu_DDE As Double
    
    Al_Cu = Al / Cu
    
c0 = 0.000285123
c1 = 1.91034
    
    photon_Cu_DDE = ((Abs((c0 + (c1 * ((dblAlCu) ^ 2))))) ^ (0.5)) * Cu
    
    jboyd_calc_photon_Cu_DDE = photon_Cu_DDE

End Function

Public Function jboyd_calc_photon_PL_DDE(PL As Double, Cu As Double)

    Dim PL_Cu As Double
    Dim photon_PL_DDE As Double

    PL_Cu = PL / Cu
    
    Dim c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, pct As Double
    pct = 1#
c0 = 1.30858
c1 = 2.96486
c2 = 0.460085
c3 = -0.0140062
c4 = 0.000121276
c5 = -0.0181811
c6 = 1.4691
c7 = -0.0411886
c8 = 0.000243623
c9 = 0.00000347525
    
    photon_PL_DDE = ((c0 + _
                (c1 * (PL_Cu)) + _
                (c2 * ((PL_Cu) ^ 2)) + _
                (c3 * ((PL_Cu) ^ 3)) + _
                (c4 * ((PL_Cu) ^ 4))) / _
                (1# + (c5 * (PL_Cu)) + _
                (c6 * ((PL_Cu) ^ 2)) + _
                (c7 * ((PL_Cu) ^ 3)) + _
                (c8 * ((PL_Cu) ^ 4)) + _
                (c9 * ((PL_Cu) ^ 5)))) * PL
                
                
    jboyd_calc_photon_PL_DDE = photon_PL_DDE
End Function

Public Function jboyd_calc_average_photon_DDE(OW As Double, PL As Double, Al As Double, Cu As Double, Optional input_photon_rad_qual As String = "all", _
            Optional input_beta_ind As Variant = -99, Optional input_beta_is_mixed As Variant = -99)
    'Options for input_photon_rad_qual are:
    ' all : coefficients calibrated from datasets considering ALL energies, i.e., PL, PM, PH
    ' PL : coefficients calibrated from datasets only with PL
    ' PM ""
    ' PH ""
    Dim photon_Cu, photon_PL As Double
    Dim average_photon_DDE As Double
    Dim beta_SDE As Double
    Dim beta_is_mixed As Integer
    
    If (input_beta_ind = -99) Then
        beta_ind = jboyd_source_type_indicator(PL, Cu, "beta_indicator")
    Else
        beta_ind = CInt(input_beta_ind)
    End If 

    If (input_beta_is_mixed = -99) Then
        beta_is_mixed = jboyd_source_type_indicator(PL, Cu, "mixed")
    Else 
        beta_is_mixed = CInt(input_beta_is_mixed)
    End If 
    
    If (StrComp("PL", input_photon_rad_qual, vbBinaryCompare) = 0) Then
        photon_PL = jboyd_calc_photon_PL_DDE_energyPL(PL, Cu)
        photon_Cu = jboyd_calc_photon_Cu_DDE_energyPL(Al, Cu)
    ElseIf (StrComp("PM", input_photon_rad_qual, vbBinaryCompare) = 0) Then
        photon_PL = jboyd_calc_photon_PL_DDE_energyPM(PL, Cu)
        photon_Cu = jboyd_calc_photon_Cu_DDE_energyPM(Al, Cu)
    ElseIf (StrComp("PH", input_photon_rad_qual, vbBinaryCompare) = 0) Then
        photon_PL = jboyd_calc_photon_PL_DDE_energyPH(PL, Cu)
        photon_Cu = jboyd_calc_photon_Cu_DDE_energyPH(Al, Cu)
    Else
        photon_PL = jboyd_calc_photon_PL_DDE(PL, Cu)
        photon_Cu = jboyd_calc_photon_Cu_DDE(Al, Cu)
    End If

    average_photon_DDE = (photon_Cu + photon_PL) / 2#
    
    If (beta_ind = 1) And (beta_is_mixed = 0) Then
        jboyd_calc_average_photon_DDE = jboyd_calc_beta_DDE(OW, PL, Al, Cu, -99, beta_is_mixed)
    Else
        jboyd_calc_average_photon_DDE = average_photon_DDE
    End If
    
    Dim c1, c2, c3, c4, calc_DDE As Double
c1 = -0.116
c2 = 0.267
c3 = 0.117
c4 = 0.646

'c1 = -0.121
'c2 = 0.275
'c3 = 0.114
'c4 = 0.658

    calc_DDE = c1 * OW + c2 * PL + c3 * Al + c4 * Cu

    jboyd_calc_average_photon_DDE = calc_DDE
    
End Function

Function jboyd_calc_SDE(OW As Double, PL As Double, Al As Double, Cu As Double, Optional input_beta_indicator As Variant = -99, _ 
            Optional input_pure_beta As Variant = -99, Optional input_beta_mixed As Variant = -99, Optional input_pure_photon As Variant = -99, _
            Optional input_gen_radiation_quality As Variant = -99, Optional input_photon_radiation_quality As Variant = -99, Optional input_high_energy_beta As Variant = -99)

    Dim total_SDE As Double
    Dim beta_SDE, photon_SDE, photon_DDE_ave, photon_DDE_rq As Double
    Dim beta_indicator, pure_beta, high_energy_beta, pure_photon As Integer
    Dim gen_radiation_quality, photon_radiation_quality As String
    
    Dim OW_PL, OW_Al, OW_Cu, PL_Cu, OWCu_PLCu_ratio, OW_Al_Cu_sub As Double
    OW_PL = OW / PL
    OW_Al = OW / Al
    OW_Cu = OW / Cu
    PL_Cu = PL / Cu
    OWCu_PLCu_ratio = (OW - Cu) / (PL - Cu)
    OW_Al_Cu_sub = (OW_Al) / (OWCu_PLCu_ratio)
    
    If (input_beta_indicator = -99) Then
        beta_indicator = jboyd_source_type_indicator(PL, Cu, "beta_indicator")
    Else
        beta_indicator = CInt(input_beta_indicator)
    End If
        
    If (input_pure_beta = -99) Then
        pure_beta = jboyd_source_type_indicator(PL, Cu, "betas")
    Else
        pure_beta = CInt(input_pure_beta)
    End If
    
    If (input_beta_mixed = -99) Then
        beta_mixed = jboyd_source_type_indicator(PL, Cu, "mixed")
    Else 
        beta_mixed = CInt(input_beta_mixed)
    End If 

    If (input_pure_photon = -99) Then 
        pure_photon = jboyd_source_type_indicator(PL, Cu, "photons")
    Else 
        pure_photon = CInt(input_pure_photon)
    End If

    If (input_gen_radiation_quality = -99) Then
        gen_radiation_quality = jboyd_calc_rad_quality(OW, PL, Al, Cu, pure_beta, beta_mixed, pure_photon)
    Else 
        gen_radiation_quality = input_gen_radiation_quality
    End If

    If (input_photon_radiation_quality = -99) Then
        photon_radiation_quality = jboyd_calc_photon_rad_quality(OW, PL, Al, Cu, 0)
    Else 
        photon_radiation_quality = input_photon_radiation_quality
    End If

    If (input_high_energy_beta = -99) Then 
        high_energy_beta = jboyd_mixed_beta_high_low_test(OW, Al) '0 if BL and 1 if BH
    Else 
        high_energy_beta = CInt(input_high_energy_beta)
    End If

    beta_SDE = jboyd_calc_beta_SDE(OW, PL, Al, Cu)
    photon_DDE_ave = jboyd_calc_average_photon_DDE(OW, PL, Al, Cu, "all", beta_indicator, beta_mixed)
    photon_DDE_rq = jboyd_calc_average_photon_DDE(OW, PL, Al, Cu, photon_radiation_quality, beta_indicator, beta_mixed)
    photon_SDE = jboyd_calc_average_photon_SDE(OW, PL, Al, Cu, gen_radiation_quality, photon_radiation_quality, pure_beta, beta_mixed, pure_photon, high_energy_beta)

    If (photon_SDE < 0) Then
        jboyd_calc_SDE = photon_DDE_rq
    Else
        jboyd_calc_SDE = photon_SDE
    End If

End Function

Public Function jboyd_calc_photon_rad_quality(OW As Double, PL As Double, Al As Double, Cu As Double, Optional return_type As Integer = 1, Optional rad_qual_mixed_input As Double = -99)
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
        PL = PL * 1.00000001
    Else
        OW_PL_cuCorr = (OW - Cu) / (PL - Cu)
    End If
    
    'Calculate Effective Energy
    If (PL_Cu) <= 1.05 Then
        energy = 662#
    Else
        energy = jboyd_calc_effective_energy(PL, Cu)
    End If
    

    If (energy >= 662) Then
        photon_source = 1
        photon_source_string = "PH"
    ElseIf (energy <= 40) Then
        photon_source = 3
        photon_source_string = "PL"
    ElseIf (40 < energy <= 200) Then
        photon_source = 2
        photon_source_string = "PM"
    Else
        photon_source = 1
        photon_source_string = "PH"
    End If
      
    If (return_type = 1) Then
        'calc_radiation_quality = beta_low_high
        jboyd_calc_photon_rad_quality = photon_source
    Else
        'calc_radiation_quality = beta_low_high_string
        jboyd_calc_photon_rad_quality = photon_source_string
    End If
    
End Function



Function jboyd_calc_luxel(OW As Double, PL As Double, Al As Double, Cu As Double, Optional return_value As String = "SDE")
    'return_value let's you select which variable to return from this function
    ' "SDE" is default and is the TOTAL SDE (photon SDE + beta SDE)
    ' "beta_SDE" just the beta_SDE value
    ' "photon_SDE" just the photon_SDE value
    ' "photon_DDE is the photon_DDE value
    ' "DDE" is total DDE
    ' "radiation_quality" is the radiation quality
    ' "beta_indicator" returns the beta indicator (after process and checks)
    ' "mixed_beta" is the boolean for whether the source is pure photon or mixed photons/betas

    Dim luxel_calc As Variant
    Dim SDE_val, beta_SDE, photon_SDE As Double
    Dim SDE_val_init As Double
    Dim DDE_val, beta_DDE, photon_DDE As Double
    Dim beta_indicator, mixed_beta, pure_beta, lin_beta_test, beta_low_high As Integer
    Dim rad_quality As String
    
    'Characterize the source and determine if we have pure beta, mixed, or pure photons.
    beta_indicator = jboyd_source_type_indicator(PL, Cu, "beta_indicator")
    pure_beta = jboyd_source_type_indicator(PL, Cu, "betas")
    mixed_beta = jboyd_source_type_indicator(PL, Cu, "mixed")
    pure_photon = jboyd_source_type_indicator(PL, Cu, "photons")
    beta_low_high = jboyd_calc_beta_rad_quality(OW, PL, Al, Cu, 1)
    
    'Radiation Quality'
    rad_quality = jboyd_calc_rad_quality(OW, PL, Al, Cu, pure_beta, mixed_beta, pure_photon)
    photon_rad_quality = jboyd_calc_photon_rad_quality(OW, PL, Al, Cu, 0)

    If (StrComp("SDE", return_value, vbBinaryCompare) = 0) Or _
        (StrComp("total_SDE", return_value, vbBinaryCompare) = 0) Or _
        (StrComp("combined_SDE", return_value, vbBinaryCompare) = 0) Then
        SDE_val = jboyd_calc_combined_SDE(OW, PL, Al, Cu)
        jboyd_calc_luxel = SDE_val
        
    ElseIf (StrComp("beta_SDE", return_value, vbBinaryCompare) = 0) Then
        If (beta_indicator = 0) Then
            beta_SDE = 0
        Else
            beta_SDE = jboyd_calc_beta_SDE(OW, PL, Al, Cu)
        End If
        jboyd_calc_luxel = beta_SDE

    ElseIf (StrComp("beta_DDE", return_value, vbBinaryCompare) = 0) Then
        beta_DDE = jboyd_calc_beta_DDE(OW, PL, Al, Cu, beta_low_high, mixed_beta)
        jboyd_calc_luxel = beta_DDE

    ElseIf (StrComp("photon_SDE", return_value, vbBinaryCompare) = 0) Then
        photon_SDE = jboyd_calc_SDE(OW, PL, Al, Cu, beta_indicator, pure_beta, mixed_beta, pure_photon, rad_quality, photon_rad_quality)
        jboyd_calc_luxel = photon_SDE
               
    ElseIf (StrComp("photon_DDE", return_value, vbBinaryCompare) = 0) Then
        DDE_val = jboyd_calc_average_photon_DDE(OW, PL, Al, Cu, rad_quality, beta_indicator, mixed_beta)
        jboyd_calc_luxel = DDE_val
        
    ElseIf (StrComp("DDE", return_value, vbBinaryCompare) = 0) Or _ 
            (StrComp("total_DDE", return_value, vbBinaryCompare) = 0) Then
        If (pure_beta = 1) Then
            DDE_val = jboyd_calc_beta_DDE(OW, PL, Al, Cu, beta_low_high, mixed_beta)
        ElseIf (mixed_beta = 1) Then
            DDE_val = jboyd_calc_mixed_DDE(OW, PL, Al, Cu)
        Else 'Pure Photons'
            DDE_val = jboyd_calc_average_photon_DDE(OW, PL, Al, Cu, rad_quality, beta_indicator, mixed_beta)
        End If
        
        jboyd_calc_luxel = DDE_val
        
    ElseIf (StrComp("radiation_quality", return_value, vbBinaryCompare) = 0) Then
        jboyd_calc_luxel = rad_quality
    ElseIf (StrComp("beta_indicator", return_value, vbBinaryCompare) = 0) Then
        jboyd_calc_luxel = beta_indicator
    ElseIf (StrComp("mixed_beta", return_value, vbBinaryCompare) = 0) Then
        jboyd_calc_luxel = mixed_beta
    Else
        SDE_val = jboyd_calc_combined_SDE(OW, PL, Al, Cu)
        jboyd_calc_luxel = SDE_val
    End If
        
End Function





