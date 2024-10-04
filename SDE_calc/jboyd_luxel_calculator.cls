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
    
    ' c1 = 0.406
    ' c2 = 0.764
    ' c3 = -1.47
    ' c4 = 1.285

'Updated with new data (09/09/2024)'
' c1 = 0.26421
' c2 = 0.55419
' c3 = -0.73292
' c4 = 0.87952

'Updated against correct Req_SDE values: 17/09/2024'
' c1 = 0.27108
' c2 = 0.11047
' c3 = -0.27306
' c4 = 0.68858

'Updated with cleaned dataset from M. Salanksy, 01/10/2024'
c1 = 0.30606
c2 = 0.11738
c3 = -0.30321
c4 = 0.76364

'12/09/2024: Coefficients with Angled Data'
' c1 = 0.02311
' c2 = 0.35837
' c3 = -0.09815
' c4 = 0.46013
    combined_SDE = c1 * OW + c2 * PL + c3 * Al + c4 * Cu

    jboyd_calc_photons_combined_SDE = combined_SDE
End Function

' 'Out of date - calibrated with OLD data'
' Function jboyd_calc_mixed_combined_SDE(OW As Double, PL As Double, Al As Double, Cu As Double)

'     'mixed_beta_indicator = 0 --> NOT MIXED
'     'mixed_beta_indicator = 1 --> MIXED
'     Dim mixed_beta_indicator As Integer
'     Dim c0, c1, c2, c3 As Double
    
'     Dim OW_Cu, PL_Cu As Double

'     If (PL = Cu) Then
'         PL = 1.000001 * PL
'     End If
'     PL_Cu = PL / Cu
'     OW_Cu = OW / Cu
    
'     'I think we can get a more accurate calculation if we can distinguish between M30 and the others.
'     'Let's try to isolate M30 by checking OW/Cu andPL/Cu
'     'We will use linear functions from SDE vs (OW + PL + Al + Cu):
'     Dim m30_m, m30_b, Sum As Double
'     Dim all_m, all_b, combined_SDE, all_sum As Double
    
' 'Standard Coefficients Method:

'     If (OW_Cu > 22) And (PL_Cu > 22) Then  'M30
        
'         'Coefficient Method:
'         c0 = 0.43101
'         c1 = 0.0619596
'         c2 = -0.533692
'         c3 = 9.04421

'         combined_SDE = 1.08 * (c0 * OW + c1 * PL + c2 * Al + c3 * Cu)
' '
'         'Line method:
' '        m30_m = 0.269851
' '        m30_b = 321.049
' '        sum = OW + PL + Al + Cu
' '
' '        combined_SDE = m30_m * sum + m30_b
        
'     Else 'All Others
'         c0 = 1.48096
'         c1 = -0.878005
'         c2 = -0.716098
'         c3 = 1.23318
        
' '        combined_SDE = 1.25 * (c0 * OW + c1 * PL + c2 * Al + c3 * Cu)
        
'         combined_SDE = (OW + PL + Al + Cu) * 0.0728819 + 404.345

'         'Line method:
' '        m30_m = 0.117222
' '        m30_b = 362.638
' '        sum = OW + PL + Al + Cu
' '
' '        combined_SDE = m30_m * sum + m30_b
'     End If
    
'     jboyd_calc_mixed_combined_SDE = combined_SDE
' End Function

Function jboyd_calc_mixed_combined_SDE(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim OW_Cu, PL_Cu As Double
    Dim c1, c2, c3, c4 As Double
    Dim multiplier As Double

    multiplier = 1.0

    If (PL = Cu) Then
        PL = 1.000001*PL
    End If 

    OW_Cu = OW/Cu 
    PL_Cu = PL/Cu 

    'Different coefficients for each energy range (PL, PM, PH)'

'PH'
    If (PL_Cu <= 1.05) Then
        ' c1 = 1.41095
        ' c2 = -1.25186
        ' c3 = 0.05728
        ' c4 = 0.72176
        ' multiplier = 1.15       

    'Updated with *cleaned* dataset from M. Salanksy, 01/10/2024'
        c1 = 2.97673
        c2 = -1.04373
        c3 = -0.24034
        c4 = -0.63470
        multiplier = 1.15
    Else

'PL and PM'
'Here we can't really identify between Mixed PL and PM
'So, we can only make exceptions on particular sources like NS20, L50, M30, NS30
        If (OW_Cu > 60) And (PL_Cu > 50) Then 'Mixed NS20 (BL) Exception'
            'Updated with proper Req_SDE values: 17/09/2024'
            'Calibrated with only mixed NS20'
            c1 = 0.81314
            c2 = -0.65292
            c3 = -0.00834
            c4 = 8.32250
            multiplier = 1.06

        ' ElseIf (OW_Cu >= 22) And (OW_Cu <= 38.6) And (PL_Cu >= 22) And (PL_Cu <= 33.7) Then
        ElseIf (OW_Cu >= 16) And (OW_Cu <= 38.6) And (PL_Cu >= 16) And (PL_Cu <= 33.7) Then

            'L50 and NS30 and M50 Exception'
            'These are Low Energy Photons (PL) and not used in the ANSI Mixture consideration'
            'Calibrated with only 0deg Mixed L50 and NS30'
            ' c1 = 1.69282
            ' c2 = -0.99456
            ' c3 = -0.31328
            ' c4 = -6.31812
            'added m50'
            c1 = 1.24480
            c2 = -0.87158
            c3 = -0.10547
            c4 = -1.07830
            multiplier = 0.95

        ElseIf (OW_Cu >= 31) And (OW_Cu <= 63) And (PL_Cu >= 27) And (PL_Cu <= 56) Then
            'M30 Exception'
            'These are Low Energy Photons (PL) and not used in the ANSI Mixture consideration'
            'Calibrated with only 0deg Mixed M30'
            c1 = 0.93589
            c2 = -0.46445
            c3 = -0.47781
            c4 = 0.73547
            multiplier = 1.00

        Else 
            'This is for everything else'
            'Calibrated using Medium and High Energy (>70 keV) Angled Data
            '0deg Mixed data (except NS20, L50, M30, NS30)' 
            ' c1 = 1.57765
            ' c2 = -1.11465
            ' c3 = -0.36927
            ' c4 = 0.89249
'Calibrated without M50':
            c1 = 1.60720
            c2 = -1.14241
            c3 = -0.36884
            c4 = 0.88533

            multiplier = 1.03
        End If 

'Updated with *cleaned* dataset from M. Salasky, 02/10/2024'
'Attempting single function for all mixed sources'
c1 = 2.22057
c2 = -2.11817
c3 = 0.08370
c4 = 1.00112
multiplier = 1.0

    End If

' c1 = 2.69527
' c2 = -2.69090
' c3 = 0.27367
' c4 = 0.86661
' multiplier = 1.0

    jboyd_calc_mixed_combined_SDE = multiplier*(c1*OW + c2*PL + c3*Al + c4*Cu)

End Function

' 'OutDated -- Used old data and old methods'
' Public Function jboyd_calc_beta_rad_quality(OW As Double, PL As Double, Al As Double, Cu As Double, Optional return_type As Integer = 1, Optional input_pure_beta As Variant = -99, Optional rad_qual_mixed_input As Variant = -99)
'     'Trying to differentiate between Sr and Kr-85 here.
'     'Kr-85 is lower energy but has larger ratio values
    
'     'Beta Low = 0
'     'Beta High = 1
    
'     'return_type = 1 --> Return an integer
'     'return_type = 0 --> Return String
    
'     Dim beta_low_high As Integer
'     beta_low_high = 0
'     Dim beta_low_high_string As String
'     beta_low_high_string = ""
    
'     Dim beta_indicator As Integer
'     Dim pure_beta As Integer
    
'     If (input_pure_beta = -99) Then
'         pure_beta = jboyd_source_type_indicator(PL, Cu, "betas")
'     Else 
'         pure_beta = CInt(input_pure_beta)
'     End If
    
'     'Mixed sources
'     Dim mixed_beta As Integer
'     'If no input (mixed_beta_input) is provided, then we look up the mixed_beta condition.
'     'mixed_beta = 0 --> NOT mixed
'     'mixed_beta = 1 --> Mixed
'     If (rad_qual_mixed_input = -99) Then
' '        mixed_beta = jboyd_beta_indicator(OW, PL, Al, Cu, 1)
'         mixed_beta = jboyd_source_type_indicator(PL, Cu, "mixed")
'     Else
'         mixed_beta = CInt(rad_qual_mixed_input)
'     End If
    
'     Dim signal_sum As Double
    
'     Dim bh_sde, bl_sde As Double
    
'     Dim OW_PL As Double
'     Dim OW_Al As Double
'     Dim OW_Cu As Double
'     Dim PL_Cu As Double
'     Dim OW_PL_cuCorr As Double
    
'     OW_PL = OW / PL
'     OW_Al = OW / Al
'     OW_Cu = OW / Cu
'     PL_Cu = PL / Cu
    
'     If (PL = Cu) Then
'         PL = PL * 1.00000001
'     Else
'         OW_PL_cuCorr = (OW - Cu) / (PL - Cu)
'     End If
    
'     signal_sum = OW + PL + Al + Cu
    
'     Dim OW_PL_SumScaled As Double
'     Dim OW_Al_SumScaled As Double
'     Dim OW_Cu_SumScaled As Double
'     Dim PL_Cu_SumScaled As Double
'     Dim OW_PL_CuCorr_SumScaled As Double
    
'     OW_PL_SumScaled = 1000 * OW_PL / signal_sum
'     OW_Al_SumScaled = 1000 * OW_Al / signal_sum
'     OW_Cu_SumScaled = 1000 * OW_Cu / signal_sum
'     PL_Cu_SumScaled = 1000 * PL_Cu / signal_sum
'     OW_PL_CuCorr_SumScaled = 1000 * OW_PL_cuCorr / signal_sum
    
'     Dim OW_PL_max As Double
'     Dim OW_Al_max As Double
'     Dim OW_Cu_max As Double
'     Dim PL_Cu_max As Double
'     Dim OW_PL_CuCorr_max As Double
    
' '****************
' 'NOTE, WE CAN'T CLEANLY DISTINGUISH ALL LOW ENERGY BETA AND HIGH ENERGY BETA IN MIXED SOURCES
' 'There does seem to be a boundary between the two though in the OW element.
' 'We can draw a linear boundary between PL/Cu and OW and
' 'We can use a relation similar to: y = 0.057x - 28.5 (where y = PL/Cu and x = OW)
' 'If PL/Cu is larger than this value it will lie to the "left" of this line and most likely be a low-energy beta.
' 'Low energy betas will also exist above this but as of yet, there is no way to differentiate further

' ' BM = Beta Mixed (GENERAL)
' ' BLM = Low Energy Beta, Mixed
' ' BHM = High Energy Beta, Mixed

'     Dim PL_Cu_linear_test As Double
'     PL_Cu_linear_test = 0.057 * OW - 28.5
    
'     If (OW_PL_SumScaled > 5) And (OW_Al_SumScaled > 5) And (OW_Cu_SumScaled > 5) Then
'         beta_low_high = 0
'             beta_low_high_string = "BL"
'     ElseIf ((PL_Cu >= 45) And (OW_Al_SumScaled > 0.2) And (OW_Cu > 45)) Or (OW_PL_CuCorr_SumScaled <= 0) Or (OW_PL > 20) Then
'         beta_low_high = 0
'         If (mixed_beta) And (pure_beta = 0) Then
'             beta_low_high_string = "BLM"
' '        ElseIf (OW_PL > 3) And (OW_Al > 3) And (OW_Cu > 3) Then

'         Else
'             beta_low_high_string = "BL"
'         End If
'     Else
'         beta_low_high = 1
'         If (mixed_beta) And (pure_beta = 0) Or ((PL_Cu >= 45) And (OW_Al_SumScaled > 0.2) And (OW_Cu > 35)) Or (OW_PL_CuCorr_SumScaled <= 0) Or (OW_PL > 20) And (pure_beta = 0) Then
'             beta_low_high_string = "BHM"
'         ElseIf (OW_Cu > 3) And (PL_Cu > 2) Then
'             beta_low_high_string = "BH"
'         Else
'             beta_low_high_string = "Bmix"
'         End If
'     End If


'     If (return_type) Then
'         jboyd_calc_beta_rad_quality = beta_low_high
'     Else
'         jboyd_calc_beta_rad_quality = beta_low_high_string
'     End If

' End Function

Public Function jboyd_calc_beta_rad_quality(OW As Double, PL As Double, Al As Double, Cu As Double, Optional return_type As Integer = 1)
    'Trying to differentiate between Sr and Kr-85 here.
    'Kr-85 is lower energy but has larger ratio values
    'Based on "Source Indicator Plot" of PL vs Cu we can see distinct regions for BL and BH'    

    'We take the max PL value for our Kr85 dataset (16.79) and the min PL value for the Sr90 data (72.2)'
    'The Average value of this (44.5) will be our selection boundary'

    If (PL < 44.5) Then 
        'BL'
        beta_low_high = 0
        beta_low_high_string = "BL"
    Else
        'Case for PL > 166'
        'BH'
        beta_low_high = 1
        beta_low_high_string = "BH"
    End If 

    If (return_type = 1) Then
        jboyd_calc_beta_rad_quality = beta_low_high
    Else
        jboyd_calc_beta_rad_quality = beta_low_high_string
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
        rad_qual = jboyd_calc_beta_rad_quality(OW, PL, Al, Cu, 0)
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

Public Function jboyd_calc_general_DDE(OW As Double, PL As Double, Al As Double, Cu As Double)
    'Updated 09/09/2024'
    'Calibrated off of all Beta data (BL & BH)'
    Dim DDE As Double
    Dim c1 As Double
    Dim c2 As Double
    Dim c3 As Double
    Dim c4 As Double

    ' c1 = -0.05599
    ' c2 = 0.16381
    ' c3 = 0.14549
    ' c4 = 0.63064

'Updated with Angled Data:'
' c1 = 0.06422
' c2 = -0.03021
' c3 = 0.20945
' c4 = 0.62229

'New Calibrations 17/09/2024'
' c1 = -0.08091
' c2 = 0.14746
' c3 = 0.08176
' c4 = 0.67505

'Updated with *cleaned* dataset from M. Salanksy on 02/10/2024'
'Calibrated using All Pure Photons and All Mixture types (SDE and DDE)'
' c1 = -0.06115
' c2 = 0.13100
' c3 = 0.08992
' c4 = 0.73281

c1 = -0.12025
c2 = 0.20098
c3 = 0.07713
c4 = 0.73723

    DDE = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)

    jboyd_calc_general_DDE = DDE
End Function

Public Function jboyd_calc_general_SDE(OW As Double, PL As Double, Al As Double, Cu As Double, Optional input_beta_low_high As Variant = -99, Optional input_beta_mixed As Variant = -99)
    'Updated 09/09/2024'
    'Calibrated off of all Beta data (BL & BH)'
    Dim SDE As Double
    Dim c1 As Double
    Dim c2 As Double
    Dim c3 As Double
    Dim c4 As Double

    ' c1 = 1.07430
    ' c2 = -0.40279
    ' c3 = -0.51913
    ' c4 = 0.91254

'Updated using Angled Data'
c1 = 1.14483
c2 = -0.72830
c3 = -0.00834
c4 = 0.24827

    SDE = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)

    jboyd_calc_general_SDE = SDE
End Function

Public Function jboyd_calc_general_beta_DDE(OW As Double, PL As Double, Al As Double, Cu As Double, Optional input_beta_low_high As Variant = -99, Optional input_beta_mixed As Variant = -99)
    'Updated 09/09/2024'
    'Calibrated off of all Beta data (BL & BH)'
    Dim beta_DDE As Double
    Dim c1 As Double
    Dim c2 As Double
    Dim c3 As Double
    Dim c4 As Double

    c1 = -0.00021
    c2 = 0.02850
    c3 = -0.04035
    c4 = 1.01335

    beta_DDE = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)

    jboyd_calc_general_beta_DDE = beta_DDE
End Function

' Public Function jboyd_calc_general_beta_SDE(OW As Double, PL As Double, Al As Double, Cu As Double, Optional input_beta_low_high As Variant = -99, Optional input_beta_mixed As Variant = -99)
'     'Updated 09/09/2024'
'     'Calibrated off of all Beta data (BL & BH)'
'     Dim beta_SDE As Double
'     Dim c1 As Double
'     Dim c2 As Double
'     Dim c3 As Double
'     Dim c4 As Double

'     c1 = 0.82413
'     c2 = -0.90493
'     c3 = -0.15543
'     c4 = 1.61282

'     beta_SDE = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)

'     jboyd_calc_general_beta_SDE = beta_SDE
' End Function

Public Function jboyd_calc_beta_DDE(OW As Double, PL As Double, Al As Double, Cu As Double, Optional input_beta_low_high As Variant = -99)
    'beta_low_high = 99 -> Nothing specificed, will look it up
    'beta_low_high = 0 --> low energy beta
    'beta_low_high = 1 --> high energy beta
    Dim beta_low_high, beta_mixed As Integer

    If (beta_low_high = -99) Then
        beta_low_high = jboyd_calc_beta_rad_quality(OW, PL, Al, Cu, 1)
    Else 
        beta_low_high = CInt(input_beta_low_high)
    End If
    
'    Can bake in a check that has the algorithm select the highest value of SDE for worst-case scenarious
    Dim BH_DDE As Double
    Dim BL_DDE As Double
    Dim BMgen_DDE As Double
    
    Dim beta_DDE As Double
    Dim c1 As Double
    Dim c2 As Double
    Dim c3 As Double
    Dim c4 As Double

    Dim pct As Double
    pct = 1#

    BL_DDE = 0
    BH_DDE = 0

    If (beta_low_high = 0) Then
        'low energy betas
        ' c1 = 0.000285893
        ' c2 = 0.00953786
        ' c3 = 0.00899577
        ' c4 = 0.979589
'Updated 09/09/2024'
        c1 = -0.00060
        c2 = 0.02781
        c3 = 0.03093
        c4 = 0.97088
        BL_DDE = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)
    Else
        'Should be for (beta_low_high = 1)'
        'high energy betas
        ' c1 = -0.217875
        ' c2 = 0.361691
        ' c3 = 0.0944204
        ' c4 = 0.633397
'Updated 09/09/2024'
        c1 = -0.05867
        c2 = 0.11408
        c3 = -0.00636
        c4 = 0.93872

        BH_DDE = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)
    End If

    'Condition to check for pure beta --> Hard to identify pure from mixtures so we can look at total DDE calc and make some decisions
    'If the energies are 0 or less than 0 the we don't want that. We will just default it to ta standard low energy beta DDE val:
    
    jboyd_calc_beta_DDE = Largest(BL_DDE, BH_DDE)
    
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

Function jboyd_calc_general_beta_SDE(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c1, c2, c3, c4 As Double
    Dim Sum As Double
    Sum = OW + PL + Al + Cu
    'Calibrated off of mixed sources only to get a Generalized fit for Beta SDE'
'Updated 09/09/2024'
c1 = 1.09199
c2 = -1.28670
c3 = 0.28994
c4 = 0.06733

        general_beta_SDE = (c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu)
        
        jboyd_calc_general_beta_SDE = general_beta_SDE

End Function

Function jboyd_calc_BL_SDE(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c1, c2, c3, c4 As Double
    Dim Sum As Double
    Sum = OW + PL + Al + Cu
    
    Dim pure_BL, angled_BL As Boolean
    Dim BL_SDE_multiplier, BH_SDE_Multiplier As Double
    pure_BL = false 
    angled_BL = false

    'Can try and distinguish between pure betas and angled betas'
    If (PL < 44.5 And Cu <= 10) Then
        pure_BL = True
        BL_SDE_multiplier = 1.0
    End If 

    If (PL < 44.5 And Cu > 10) Then 
        angled_BL = True 
        BL_SDE_multiplier = 1.0
    End If 
'Updated with Angled Data'
' c1 = 1.85026
' c2 = 7.89371
' c3 = 5.09608
' c4 = 6.81685


    'low energy betas
        'pure
       
'BASIC FUNCTION
'        c1 = 2.28162
'        c2 = -2.06136
'        c3 = 0#
'        c4 = 0.0116431
' c1 = 1.535
' c2 = 4.628
' c3 = 25.375
' c4 = -22.764

'Updated 09/09/2024'
' c1 = 2.29523
' c2 = -2.21141
' c3 = 0.00689
' c4 = 0.09998

'Updated with Angled Data'
' c1 = 1.02315
' c2 = -0.51987
' c3 = -0.33103
' c4 = 0.81416

'Updated with *clean* dataset from M. Salasky, 02/10/2024'
c1 = 2.64734
c2 = 3.83965
c3 = 20.27848
c4 = 1.97002

        bl_sde = BL_SDE_multiplier*((c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu))
        
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

Function jboyd_calc_BH_SDE(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c1, c2, c3, c4 As Double
    Dim Sum As Double
    Sum = OW + PL + Al + Cu

    Dim pure_BH, angled_BH As Boolean
    pure_BH = false 
    angled_BH = false 
    'Can try and distinguish between pure betas and angled betas'   
    If (PL >= 44.5 And PL < 280) Then 
        angled_BH = True 
        BH_SDE_Multiplier = 1.0
    End If 

    If (PL >= 280) Then 
        pure_BH = True 
        BH_SDE_Multiplier = 1.0
    End If

'Updated with Angled Data'
' c1 = 0.59450
' c2 = -0.58955
' c3 = -0.88556
' c4 = 19.97215


    'high energy betas
        'pure

'Updated 09/09/2024'
'BASIC FUNCTION
    ' c1 = 1.51442
    ' c2 = -1.63437
    ' c3 = 0.38618
    ' c4 = -0.29304

'Updated on 02/10/2024 with *clean* dataset from M. Salasky'
    c1 = 0.19866
    c2 = 0.94002
    c3 = 0.59200
    c4 = -0.54506

'Updated with Angled Data'
' c1 = 0.59450
' c2 = -0.58955
' c3 = -0.88556
' c4 = 19.97215

       bh_sde = BH_SDE_Multiplier*((c1 * OW) + (c2 * PL) + (c3 * Al) + (c4 * Cu))

''Irreducial Form Factor Pol3 (high energy beta - Sy ONLY)
'    c0 = -649.878
'    c1 = 0.0000000180641
'    c2 = -0.000135992
'    c3 = 0.40616
'    bh_sde = (sum - c0) * (c1 * sum * sum + c2 * sum + c3)

'Before 09/09/2024 Update, this was the method:
'ln+expo (high energy beta - Sy ONLY)
    ' c0 = 27.9222
    ' c1 = 0.0001
    ' c2 = 1.5
    ' c3 = -2500
    ' c4 = -664.077

    ' bh_sde = c0 * Math.Log(Sum / c1) + c2 * Exp((-0.5) * ((Sum - c3) / c4))

    jboyd_calc_BH_SDE = bh_sde

End Function

Function jboyd_calc_BM_SDE(OW As Double, PL As Double, Al As Double, Cu As Double)
    Dim c1, c2, c3, c4 As Double
    Dim Sum As Double
    Sum = OW + PL + Al + Cu
End Function

' Function jboyd_pure_beta_test(PL As Double, Cu As Double)
    ' OUT-DATED: BASED ON OLD PETG 0.04 DATA'
'     Dim m1, m2, b1, b2, upper_test, lower_test As Double
'     Dim pure_beta As Integer
    
'     'Based on a plot of PL vs Cu for photons, betas, and mixed, we can find some boundaries between each of these types
'     'We can isolate most of the pure betas between two lines. Values between the lines are pure beta.
'     'upper:
'     m1 = 9.2329
'     b1 = 17.5
'     upper_test = m1 * Cu + b1
    
'     'lower:
'     m2 = 9.2329
'     b2 = -250
'     lower_test = m2 * 9.2329 + b2
    
'     'I will also put a cap on Cu just so that we "truncate" the region between the two lines
'     'If it goes to infinity it may overlap with the region defined by the pol4 for pure photons
    
'     If (PL < upper_test) And (PL > lower_test) And (Cu < 200) Then
'         pure_beta = 1
'     Else
'         pure_beta = 0
'     End If
    
'     jboyd_pure_beta_test = pure_beta
    
' End Function

Function jboyd_pure_beta_test(PL As Double, Cu As Double)
    'Updated 09/09/2024'
    Dim PL_limit, PL_test, Cu_limit, Cu_test As Double
    Dim pure_beta As Integer
    
    'Based on a plot of PL vs Cu for photons, betas, and mixed, we can find some boundaries between each of these types
    'We can isolate most of the pure betas between two lines. 
    'We can create a boxed region and everything within (Or below) these limits are Pure Betas'
    'These limits were determined empirically by looking at the plot and data table'
    
    ' PL_limit = 500.0
    ' Cu_limit = 52.5

    'Updated 02/10/2024 after M. Salasky provided cleaned dataset:'
    PL_limit = 338.81
    Cu_limit = 37.895  
    
    If (PL < PL_limit And Cu < Cu_limit) Then
        pure_beta = 1
    Else
        pure_beta = 0
    End If
    
    jboyd_pure_beta_test = pure_beta
    
End Function


' Function jboyd_pure_photon_test(PL As Double, Cu As Double)
    ' OUT-DATED: BASED ON OLD PETG 0.04 DATA'
'     Dim c0, c1, c2, c3, pol4_test As Double
'     Dim pure_photon As Integer

'     'Based on a plot of PL vs Cu for photons, betas, and mixed, we can find some boundaries between each of these types
'     'The boundary between pure and (mixed or beta) can be defined by a 4th order polynomial.
'     '->Values below the pol4 are mixed or beta, and values above the pol4 are pure photons
'     ' y = A + B(x-D) + C(x-D)^2, where y = PL and x = Cu
'     'Pol4:
'     c0 = 311.6
'     c1 = -2.6
'     c2 = 0.0183
'     c3 = 534
    
'     pol4_test = c0 + c1 * (Cu - c3) + c2 * ((Cu - c3) ^ (2))
    
'     If (PL > pol4_test) Then
'         pure_photon = 1
'     Else
'         pure_photon = 0
'     End If
    
'     jboyd_pure_photon_test = pure_photon
' End Function

'Updated 09/09/2024'
Function jboyd_pure_photon_test(PL As Double, Cu As Double)
    Dim c0, c1, c2, c3, c4, c5, pol5_test As Double
    Dim pure_photon As Integer

    'Based on a plot of PL vs Cu for photons, betas, and mixed, we can find some boundaries between each of these types
    'The boundary between pure and (mixed or beta) can be defined by a 5th order polynomial and horizontal line.
    'This is a piece-wise function that has the pol5 for Cu < 512.5 and a horizontal line for Cu >= 512.5'
    '->Values below the pol5 are mixed or beta, and values above the pol4 are pure photons
    ' y = c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4 + c5*x^5, where y = PL and x = Cu
    'Pol5:
    ' c0 = 11521.3
    ' c1 = -156
    ' c2 = 1.0012
    ' c3 = -0.003141
    ' c4 = 0.00000485109
    ' c5 = -0.0000000029946

    'Updated fit from cleaned/updated data from M. Salasky on 01/10/2024'
    c0 = 10569
    c1 = -184.612
    c2 = 1.5375
    c3 = -0.00621976
    c4 = 1.21312e-005
    c5 = -9.1689e-009

'Previously used a pol5 and a horizontal line. Now (02/10/2024), just the pol5 is used' 
    pol5_test = c0 + (c1 * Cu) + c2*(Cu)^(2) + c3*(Cu)^(3) + c4*(Cu)^(4) + c5*(Cu)^(5)
    ' If (Cu < 512.5) Then
    '     pol5_test = c0 + (c1 * Cu) + c2*(Cu)^(2) + c3*(Cu)^(3) + c4*(Cu)^(4) + c5*(Cu)^(5)
    ' Else
    '     'Horizontal line at PL = 534'
    '     pol5_test = 534
    ' End If
        
    If (PL > pol5_test) Then
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
 
    Dim pure_BL, angled_BL As Boolean
    Dim BL_SDE_multiplier, BH_SDE_Multiplier As Double
    pure_BL = false 
    angled_BL = false

    'Can try and distinguish between pure betas and angled betas'
    If (PL < 44.5 And Cu <= 10) Then
        pure_BL = True
        BL_SDE_multiplier = 1.0
    End If 

    If (PL < 44.5 And Cu > 10) Then 
        angled_BL = True 
        BL_SDE_multiplier = 1.0
    End If 

    'Can try two methods for calculating the combined SDE for BL betas:
    'Standard Coefficient method:
    'SEEMS TO BE THE MORE ACCURATE METHOD
' c1 = 2.366
' c2 = -3.125
' c3 = -0.005
' c4 = 0.941

'Updated with new data (09/09/2024)'
' c1 = 2.29473
' c2 = -2.13331
' c3 = 0.05090
' c4 = 1.01285

'Updated with Angled Data'
' c1 = 1.85026
' c2 = 7.89371
' c3 = 5.09608
' c4 = 6.81685

'Updated using correct Req_SDE Values':
c1 = 2.52238
c2 = -2.04248
c3 = -0.47213
c4 = 1.37408

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

    Dim pure_BH, angled_BH As Boolean
    pure_BH = false 
    angled_BH = false 
    'Can try and distinguish between pure betas and angled betas'   
    If (PL >= 44.5 And PL < 280) Then 
        angled_BH = True 
        BH_SDE_Multiplier = 1.0
    End If 

    If (PL >= 280) Then 
        pure_BH = True 
        BH_SDE_Multiplier = 1.0
    End If

    'Standard Coefficient method:

    ' c1 = 1.487
    ' c2 = -1.534
    ' c3 = 0.318
    ' c4 = 0.658
 
'Updated with new data (09/09/2024)'   
' c1 = 1.45293
' c2 = -1.51168
' c3 = 0.37553
' c4 = 0.63825

'Updated with Angled Data'
' c1 = 0.59450
' c2 = -0.58955
' c3 = -0.88556
' c4 = 19.97215

'Updated using correct Req_SDE values:'
c1 = 0.22246
c2 = 0.75279
c3 = 0.54858
c4 = -0.45567
    

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
    Dim general_beta_SDE As Double
    Dim bh_sde As Double
    Dim bl_sde As Double
    
    Dim final_MAX_beta_SDE As Double

    bh_sde = jboyd_calc_BH_SDE(OW, PL, Al, Cu)
    bl_sde = jboyd_calc_BL_SDE(OW, PL, Al, Cu)
    general_beta_SDE = jboyd_calc_general_beta_SDE(OW, PL, Al, Cu)

    If (pure_photons = 1) Then
        jboyd_calc_beta_SDE = 0
    Else
        If (beta_low_high = 0) Then
            'low energy betas
    '        If (beta_mixed = 0) Then
            If (pure_beta = 1) Then
                final_MAX_beta_SDE = bl_sde
            Else
                final_MAX_beta_SDE = general_beta_SDE
            End If
        Else
            'high energy betas
    '        If (beta_mixed = 0) Then
            If (pure_beta = 1) Then
                final_MAX_beta_SDE = bh_sde
            Else
                final_MAX_beta_SDE = general_beta_SDE
            End If
        End If
    End If
    
    jboyd_calc_beta_SDE = final_MAX_beta_SDE
    
End Function

'Out of date -- calibrated with old data'
' Public Function jboyd_calc_photon_rad_quality(OW As Double, PL As Double, Al As Double, Cu As Double, Optional return_type As Integer = 1, Optional rad_qual_mixed_input As Double = -99)
'     Dim photon_source As Double
'     Dim photon_source_string As String
'     Dim energy As Double
            
'     photon_source = 0
'     photon_source_string = ""
'     'photon_source = 0 --> No photons
'     'photon_source = 1 --> Low CVs, High Energy, PH
'     'photon_source = 2 --> Mid-energy, PM
'     'photon_source = 3 --> Low-energy, PL
    
'     'Beta Low = 0
'     'Beta High = 1
    
'     'return_type = 1 --> Return an integer
'     'return_type = 0 --> Return String

'     Dim OW_PL As Double
'     Dim OW_Al As Double
'     Dim OW_Cu As Double
'     Dim PL_Cu As Double
'     Dim OW_PL_cuCorr As Double
    
'     OW_PL = OW / PL
'     OW_Al = OW / Al
'     OW_Cu = OW / Cu
'     PL_Cu = PL / Cu
    
'     If (PL = Cu) Then
'         PL = PL * 1.00000001
'     Else
'         OW_PL_cuCorr = (OW - Cu) / (PL - Cu)
'     End If
    
'     'Calculate Effective Energy
'     If (PL_Cu) <= 1.05 Then
'         energy = 662#
'     Else
'         energy = jboyd_calc_effective_energy(PL, Cu)
'     End If
    
'     'Calculate Radiation Quality'
'     If (energy >= 662) Then
'         photon_source = 1
'         photon_source_string = "PH"
'     ElseIf (energy <= 40) Then
'         photon_source = 3
'         photon_source_string = "PL"
'     ElseIf (40 < energy <= 200) Then
'         photon_source = 2
'         photon_source_string = "PM"
'     Else
'         photon_source = 1
'         photon_source_string = "PH"
'     End If
      
'     If (return_type = 1) Then
'         'calc_radiation_quality = beta_low_high
'         jboyd_calc_photon_rad_quality = photon_source
'     Else
'         'calc_radiation_quality = beta_low_high_string
'         jboyd_calc_photon_rad_quality = photon_source_string
'     End If
    
' End Function

Public Function jboyd_calc_photon_rad_quality(OW As Double, PL As Double, Al As Double, Cu As Double, Optional return_type As Integer = 1, Optional rad_qual_mixed_input As Double = -99)
    'Updated as of 09/09/2024'
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
        'PH'
        energy = 662#
    Else
        energy = jboyd_calc_effective_energy(PL, Cu)
    End If
    
    'Calculate Radiation Quality'
    'Instead of using the energy, we can make a determination on the value of PL_Cu just like we did with PH above'
    ' If (PL_Cu <= 1.05) Then
    '     photon_source = 1
    '     photon_source_string = "PH"
    ' ElseIf (PL_Cu > 1.05) And (PL_Cu < 6) Then
    '     photon_source = 2
    '     photon_source_string = "PM"
    ' Else
    '     'Here, we should be in the range of PL_Cu > 6'
    '     photon_source = 3
    '     photon_source_string = "PL"
    ' End If
    ' If( energy > 300) Then
    '     photon_source = 1
    '     photon_source_string = "PH"
    ' End If 

'Modified 18/09/2024'
    'We will still distinguish pure PH using the same energy discrimination'
    'Instead of using PL_Cu for the others, we can just use the value of Cu.'
    'We can look at a plot of PL vs Cu and cut on certain Cu values'
    'We can average values between boundaries:'



    Dim lowEph_Cu_max, midEph_Cu_min, lowEph_midEph_Cu_limit As Double
    ' lowEph_Cu_max = 299.6
    ' midEph_Cu_min = 340.7
    'Updated with *cleaned* dataset from M. Salanksy, 02/10/2024:'
    lowEph_Cu_max = 301.74
    midEph_Cu_min = 306.39

    lowEph_Cu_midEp_limit = (lowEph_Cu_max + midEph_Cu_min)/2

    If (PL_Cu <= 1.05) Then
        photon_source = 1
        photon_source_string = "PH"
    ElseIf (Cu <= lowEph_Cu_midEp_limit) Then
        photon_source = 2
        photon_source_string = "PL"
    Else
        photon_source = 3
        photon_source_string = "PM"
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
        SDE_val = jboyd_calc_combined_SDE(OW, PL, Al, Cu, beta_indicator, pure_beta, mixed_beta, pure_photon, beta_low_high)
        jboyd_calc_luxel = SDE_val
        
    ElseIf (StrComp("beta_SDE", return_value, vbBinaryCompare) = 0) Then
        If (beta_indicator = 0) Then
            beta_SDE = 0
        Else
            beta_SDE = jboyd_calc_beta_SDE(OW, PL, Al, Cu)
        End If
        jboyd_calc_luxel = beta_SDE

    ElseIf (StrComp("beta_DDE", return_value, vbBinaryCompare) = 0) Then
        If (pure_beta = 1) Then
            beta_DDE = jboyd_calc_beta_DDE(OW, PL, Al, Cu, beta_low_high)
        ElseIf (mixed_beta = 1) Then 
            beta_DDE = jboyd_calc_general_DDE(OW, PL, Al, Cu)
        Else 'Pure Photons'
            beta_DDE = 0
        End If
        jboyd_calc_luxel = beta_DDE

    ElseIf (StrComp("photon_SDE", return_value, vbBinaryCompare) = 0) Then
        photon_SDE = jboyd_calc_photons_combined_SDE(OW, PL, Al, Cu)
        jboyd_calc_luxel = photon_SDE
               
    ElseIf (StrComp("photon_DDE", return_value, vbBinaryCompare) = 0) Then
        If (pure_beta = 1) Then 
            DDE_val = 0
        Else 'Pure Photons and Mixed'
            DDE_val = jboyd_calc_general_DDE(OW, PL, Al, Cu)
        End If
        jboyd_calc_luxel = DDE_val
        
    ElseIf (StrComp("DDE", return_value, vbBinaryCompare) = 0) Or _ 
            (StrComp("total_DDE", return_value, vbBinaryCompare) = 0) Then
        If (pure_beta = 1) Then
            ' DDE_val = jboyd_calc_beta_DDE(OW, PL, Al, Cu, beta_low_high)
            'If we have a pure beta tthen we know that Beta DDE is zero, but we can set it to 0.001'
            DDE_val = 0.001
        ' ElseIf (mixed_beta = 1) Then
        '     DDE_val = jboyd_calc_mixed_DDE(OW, PL, Al, Cu)
        Else 'Pure Photons and Mixed Sources'
            DDE_val = jboyd_calc_general_DDE(OW, PL, Al, Cu)
        End If
        
        jboyd_calc_luxel = DDE_val
        
    ElseIf (StrComp("radiation_quality", return_value, vbBinaryCompare) = 0) Then
        jboyd_calc_luxel = rad_quality
    ElseIf (StrComp("beta_indicator", return_value, vbBinaryCompare) = 0) Then
        jboyd_calc_luxel = beta_indicator
    ElseIf (StrComp("mixed_beta", return_value, vbBinaryCompare) = 0) Then
        jboyd_calc_luxel = mixed_beta
    Else
        SDE_val = jboyd_calc_combined_SDE(OW, PL, Al, Cu, beta_indicator, pure_beta, mixed_beta, pure_photon, beta_low_high)
        jboyd_calc_luxel = SDE_val
    End If
        
End Function


