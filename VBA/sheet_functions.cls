Sub PopulateMixtures(startCell As Range, beta As Double, photon As Double)
    Dim i As Integer
    Dim mixture As Double

    ' Loop to populate the cells with different mixtures
    For i = 1 To 9
        mixture = (10 - i) / 10 * photon + i / 10 * beta
        startCell.offset(i - 1, 0).Value = mixture
    Next i
End Sub

Sub PopulateMixturesWithSelection()
    Dim startCell As Range
    Dim betaCell As Range
    Dim photonCell As Range
    Dim beta As Double
    Dim photon As Double
    Dim i As Integer
    Dim mixture As Double

    ' Prompt user to select the starting cell
    Set startCell = Application.InputBox("Select the starting cell:", Type:=8)
    
    ' Prompt user to select the cell containing the beta value
    Set betaCell = Application.InputBox("Select the cell containing the beta value:", Type:=8)
    
    ' Prompt user to select the cell containing the photon value
    Set photonCell = Application.InputBox("Select the cell containing the photon value:", Type:=8)
    
    ' Get the values from the selected cells
    beta = betaCell.Value
    photon = photonCell.Value

    ' Loop to populate the cells with different mixtures
    For i = 1 To 9
        mixture = (10 - i) / 10 * photon + i / 10 * beta
        startCell.offset(i - 1, 0).Value = mixture
    Next i
End Sub

Sub PopulateMixturesWithDetails()
    Dim startCell As Range
    Dim betaCell As Range
    Dim photonCell As Range
    Dim betaNameCell As Range
    Dim photonNameCell As Range
    Dim beta As Double
    Dim photon As Double
    Dim betaName As String
    Dim photonName As String
    Dim i As Integer
    Dim j As Integer
    Dim mixture As Double
    Dim ratio As String

    ' Prompt user to select the starting cell for numeric output
    Set startCell = Application.InputBox("Select the starting cell for numeric output:", Type:=8)
    
    ' Prompt user to select the cell containing the first beta value
    Set betaCell = Application.InputBox("Select the starting cell for beta values:", Type:=8)
    
    ' Prompt user to select the cell containing the first photon value
    Set photonCell = Application.InputBox("Select the starting cell for photon values:", Type:=8)
    
    ' Prompt user to select the cell containing the beta name
    Set betaNameCell = Application.InputBox("Select the cell containing the beta name:", Type:=8)
    
    ' Prompt user to select the cell containing the photon name
    Set photonNameCell = Application.InputBox("Select the cell containing the photon name:", Type:=8)
    
    ' Get the names from the selected cells
    betaName = betaNameCell.Value
    photonName = photonNameCell.Value

    ' Populate the columns with names and ratios on the left of the starting cell
    For i = 1 To 9
        ratio = i & ":" & (10 - i)
        startCell.offset(i, -3).Value = photonName
        startCell.offset(i, -2).Value = betaName
        startCell.offset(i, -1).Value = ratio
    Next i

    ' Loop through the blocks of mixtures
    For j = 0 To 4
        ' Loop to populate the cells with different mixtures
        For i = 9 To 1 Step -1
            ' Get the values for beta and photon for this row
            beta = betaCell.offset(0, j).Value
            photon = photonCell.offset(0, j).Value
            
            ' Calculate the mixture and populate the cell
            mixture = (10 - i) / 10 * photon + i / 10 * beta
            startCell.offset(10 - i, j).Value = mixture
        Next i
    Next j
End Sub

Sub PopulateMixturesForBetaAndPhotonsSingle()
    Dim startCell As Range
    Dim betaCell As Range
    Dim betaNameCell As Range
    Dim photonCells As Range
    Dim photonNames As Range
    Dim beta As Double
    Dim betaName As String
    Dim i As Integer
    Dim j As Integer
    Dim mixture As Double
    Dim ratio As String
    Dim photonValues() As Variant
    Dim photonNames() As Variant
    Dim lastRow As Long
    Dim currentRow As Long

    ' Prompt user to select the starting cell for numeric output
    Set startCell = Application.InputBox("Select the starting cell for numeric output:", Type:=8)
    
    ' Prompt user to select the cell containing the beta value
    Set betaCell = Application.InputBox("Select the cell containing the beta value:", Type:=8)
    
    ' Prompt user to select the cell containing the beta name
    Set betaNameCell = Application.InputBox("Select the cell containing the beta name:", Type:=8)
    
    ' Get the beta value and name
    beta = betaCell.Value
    betaName = betaNameCell.Value

    ' Prompt user to select the range of photon values
    Set photonCells = Application.InputBox("Select the range of cells containing the photon values:", Type:=8)
    
    ' Prompt user to select the range of photon names
    Set photonNames = Application.InputBox("Select the range of cells containing the photon names:", Type:=8)
    
    ' Get the values and names of photons
    photonValues = photonCells.Value
    photonNames = photonNames.Value
    
    ' Initial row for the first set of results
    currentRow = startCell.Row
    
    ' Loop through each photon value
    For j = 1 To UBound(photonValues, 1)
        ' Get the current photon value and name
        Dim photon As Double
        Dim photonName As String
        
        photon = photonValues(j, 1)
        photonName = photonNames(j, 1)
        
        ' Populate the columns with names and ratios for this set of photon values
        For i = 1 To 9
            ratio = (10 - i) & ":" & i
            startCell.offset(currentRow - startCell.Row + i - 1, -3).Value = photonName
            startCell.offset(currentRow - startCell.Row + i - 1, -2).Value = betaName
            startCell.offset(currentRow - startCell.Row + i - 1, -1).Value = ratio
        Next i

        ' Calculate the mixtures for the current photon value
        For i = 9 To 1 Step -1
            mixture = (10 - i) / 10 * photon + i / 10 * beta
            startCell.offset(currentRow - startCell.Row + 10 - i - 1, 0).Value = mixture
        Next i
        
        ' Move the output row down for the next photon value's calculations
        currentRow = currentRow + 10 ' 9 rows of mixtures + 1 row for spacing
    Next j
End Sub

'Sub PopulateMixturesForBetaAndPhotons()
'    Dim startCell As Range
'    Dim betaCell As Range
'    Dim betaNameCell As Range
'    Dim photonCells As Range
'    Dim photonCell As Range
'    Dim photonTypes As Range
'    Dim photonTypeCell As Range
'    Dim beta As Double
'    Dim betaName As String
'    Dim i As Integer, k As Integer
'    Dim mixture As Double
'    Dim ratio As String
'    Dim currentRow As Long
'    Dim photon As Double
'    Dim photonType As String
'    Dim colOffset As Integer
'    Dim calcColOffset As Integer
'
'    ' Prompt user to select the starting cell for numeric output
'    On Error Resume Next
'    Set startCell = Application.InputBox("Select the starting cell for numeric output:", Type:=8)
'    On Error GoTo 0
'    If startCell Is Nothing Then Exit Sub ' User canceled
'
'    ' Prompt user to select the cell containing the beta value
'    On Error Resume Next
'    Set betaCell = Application.InputBox("Select the cell containing the beta value:", Type:=8)
'    On Error GoTo 0
'    If betaCell Is Nothing Then Exit Sub ' User canceled
'
'    ' Prompt user to select the cell containing the beta name
'    On Error Resume Next
'    Set betaNameCell = Application.InputBox("Select the cell containing the beta name:", Type:=8)
'    On Error GoTo 0
'    If betaNameCell Is Nothing Then Exit Sub ' User canceled
'
'    ' Get the beta value and name
'    beta = betaCell.Value
'    betaName = betaNameCell.Value
'
'    ' Prompt user to select the range of photon values (allowing non-contiguous selection)
'    On Error Resume Next
'    Set photonCells = Application.InputBox("Select the cells containing the photon values (you can select multiple non-contiguous cells):", Type:=8)
'    On Error GoTo 0
'    If photonCells Is Nothing Then Exit Sub ' User canceled
'
'    ' Prompt user to select the range of photon types
'    On Error Resume Next
'    Set photonTypes = Application.InputBox("Select the cells containing the photon types (ensure they correspond to the photon values):", Type:=8)
'    On Error GoTo 0
'    If photonTypes Is Nothing Then Exit Sub ' User canceled
'
'    ' Initial row for the first set of results
'    currentRow = startCell.Row
'
'    ' First, print the ratios and initial mixtures for the first set of columns
'    For Each photonCell In photonCells
'        photon = photonCell.Value
'
'        ' Get the corresponding photon type
'        Set photonTypeCell = photonTypes.Cells(photonCell.Row - photonCells.Cells(1).Row + 1, photonCell.Column - photonCells.Cells(1).Column + 1)
'        photonType = photonTypeCell.Value
'
'        ' Populate the columns with types and ratios for this set of photon values (once)
'        For i = 1 To 9
'            ratio = i & ":" & (10 - i)
'            startCell.Offset(currentRow - startCell.Row + i - 1, -3).Value = photonType
'            startCell.Offset(currentRow - startCell.Row + i - 1, -2).Value = betaName
'            startCell.Offset(currentRow - startCell.Row + i - 1, -1).Value = ratio
'        Next i
'
'        ' Calculate the mixtures for the first photon value
'        For i = 9 To 1 Step -1
'            mixture = (10 - i) / 10 * photon + i / 10 * beta
'            startCell.Offset(currentRow - startCell.Row + 10 - i - 1, 0).Value = mixture
'        Next i
'
'        ' Move the output row down for the next photon value's calculations
'        currentRow = currentRow + 10 ' 9 rows of mixtures + 1 row for spacing
'    Next photonCell
'
'    ' Reset the currentRow for additional columns of mixtures
'    currentRow = startCell.Row
'
'    ' Loop through the next 3 sets of photon and beta values (columns)
'    For k = 1 To 4 ' 1 is the second column, 2 is the third, etc.
'        ' Adjust for the current column offset (starting from the second set of columns)
'        colOffset = k
'        calcColOffset = colOffset * 1 ' Adjust the number of columns to move to the right for each set
'
'        ' Get the current beta and beta name from the correct column
'        beta = betaCell.Offset(0, colOffset).Value
'        betaName = betaNameCell.Offset(0, colOffset).Value
'
'        ' Loop through each selected photon value and corresponding type
'        For Each photonCell In photonCells
'            photon = photonCell.Offset(0, colOffset).Value
'
'            ' Calculate the mixtures for the current photon value
'            For i = 9 To 1 Step -1
'                mixture = (10 - i) / 10 * photon + i / 10 * beta
'                startCell.Offset(currentRow - startCell.Row + 10 - i - 1, calcColOffset).Value = mixture
'            Next i
'
'            ' Move the output row down for the next photon value's calculations
'            currentRow = currentRow + 10 ' 9 rows of mixtures + 1 row for spacing
'        Next photonCell
'
'        ' Reset the currentRow for the next set of columns
'        currentRow = startCell.Row
'    Next k
'End Sub

Sub PopulateMixturesForBetaAndPhotons()
    Dim startCell As Range
    Dim betaCell As Range
    Dim betaNameCell As Range
    Dim photonCells As Range
    Dim photonCell As Range
    Dim photonNameColumn As Range
    Dim LDEColumn As Range
    Dim SDEColumn As Range
    Dim BetaSDEColumn As Range
    Dim beta As Double
    Dim betaName As String
    Dim i As Integer, k As Integer
    Dim mixture As Double
    Dim ratio As String
    Dim currentRow As Long
    Dim photon As Double
    Dim photonType As String
    Dim LDEPhotonValue As Double
    Dim SDEPhotonValue As Double
    Dim BetaSDEPhotonValue As Double
    Dim LDEBetaValue As Double
    Dim SDEBetaValue As Double
    Dim BetaSDEBetaValue As Double

    ' Prompt user to select the starting cell for numeric output
    On Error Resume Next
    Set startCell = Application.InputBox("Select the starting cell for numeric output:", Type:=8)
    On Error GoTo 0
    If startCell Is Nothing Then Exit Sub ' User canceled

    ' Prompt user to select the cell containing the beta value
    On Error Resume Next
    Set betaCell = Application.InputBox("Select the cell containing the beta value:", Type:=8)
    On Error GoTo 0
    If betaCell Is Nothing Then Exit Sub ' User canceled

    ' Prompt user to select the cell containing the beta name
    On Error Resume Next
    Set betaNameCell = Application.InputBox("Select the cell containing the beta name:", Type:=8)
    On Error GoTo 0
    If betaNameCell Is Nothing Then Exit Sub ' User canceled

    ' Get the beta value and name
    beta = betaCell.Value
    betaName = betaNameCell.Value

    ' Prompt user to select the range of photon values (allowing non-contiguous selection)
    On Error Resume Next
    Set photonCells = Application.InputBox("Select the cells containing the photon values (you can select multiple non-contiguous cells):", Type:=8)
    On Error GoTo 0
    If photonCells Is Nothing Then Exit Sub ' User canceled

    ' Prompt user to select the column containing the photon names
    On Error Resume Next
    Set photonNameColumn = Application.InputBox("Select the column containing the photon names:", Type:=8)
    On Error GoTo 0
    If photonNameColumn Is Nothing Then Exit Sub ' User canceled

    ' Prompt user to select the column for LDE data
    On Error Resume Next
    Set LDEColumn = Application.InputBox("Select the column containing the LDE data:", Type:=8)
    On Error GoTo 0
    If LDEColumn Is Nothing Then Exit Sub ' User canceled

    ' Prompt user to select the column for SDE data
    On Error Resume Next
    Set SDEColumn = Application.InputBox("Select the column containing the SDE data:", Type:=8)
    On Error GoTo 0
    If SDEColumn Is Nothing Then Exit Sub ' User canceled

    ' Prompt user to select the column for Beta SDE data
    On Error Resume Next
    Set BetaSDEColumn = Application.InputBox("Select the column containing the Beta SDE data:", Type:=8)
    On Error GoTo 0
    If BetaSDEColumn Is Nothing Then Exit Sub ' User canceled

    ' Determine the start row for calculations
    currentRow = startCell.Row

    ' First, print the ratios and initial mixtures for the first set of columns
    For Each photonCell In photonCells
        photon = photonCell.Value
        ' Use the row of photonCell to fetch the photon name
        photonType = photonNameColumn.Cells(photonCell.Row - photonCells.Cells(1).Row + 1, 1).Value

        ' Populate the columns with types and ratios for this set of photon values (once)
        For i = 1 To 9
            ratio = i & ":" & (10 - i)
            startCell.offset(currentRow - startCell.Row + i - 1, -3).Value = photonType
            startCell.offset(currentRow - startCell.Row + i - 1, -2).Value = betaName
            startCell.offset(currentRow - startCell.Row + i - 1, -1).Value = ratio
        Next i

        ' Calculate the mixtures for the first photon value
        For i = 9 To 1 Step -1
            mixture = (10 - i) / 10 * photon + i / 10 * beta
            startCell.offset(currentRow - startCell.Row + 10 - i - 1, 0).Value = mixture
        Next i

        ' Move the output row down for the next photon value's calculations
        currentRow = currentRow + 10 ' 9 rows of mixtures + 1 row for spacing
    Next photonCell

    ' Reset the currentRow for additional columns of mixtures
    currentRow = startCell.Row

    ' Loop through the next 4 sets of photon and beta values (columns)
    For k = 1 To 4 ' 1 is the second column, 2 is the third, etc.
        ' Adjust for the current column offset (starting from the second set of columns)
        colOffset = k
        calcColOffset = colOffset * 1 ' Adjust the number of columns to move to the right for each set

        ' Get the current beta and beta name from the correct column
        beta = betaCell.offset(0, colOffset).Value
        betaName = betaNameCell.offset(0, colOffset).Value

        ' Loop through each selected photon value and corresponding type
        For Each photonCell In photonCells
            photon = photonCell.offset(0, colOffset).Value

            ' Calculate the mixtures for the current photon value
            For i = 9 To 1 Step -1
                mixture = (10 - i) / 10 * photon + i / 10 * beta
                startCell.offset(currentRow - startCell.Row + 10 - i - 1, calcColOffset).Value = mixture
            Next i

            ' Move the output row down for the next photon value's calculations
            currentRow = currentRow + 10 ' 9 rows of mixtures + 1 row for spacing
        Next photonCell

        ' Reset the currentRow for the next set of columns
        currentRow = startCell.Row
    Next k

    ' Calculate LDE mixtures
    currentRow = startCell.Row
    For Each photonCell In photonCells
        photon = photonCell.Value
        ' Use the row of photonCell to fetch the LDE photon value
        LDEPhotonValue = LDEColumn.Cells(photonCell.Row - photonCells.Cells(1).Row + 1, 1).Value
        ' Fetch the LDE beta value
        LDEBetaValue = betaCell.offset(0, LDEColumn.Column - startCell.Column).Value

        ' Calculate the mixtures for the LDE column
        For i = 9 To 1 Step -1
            mixture = (10 - i) / 10 * LDEPhotonValue + i / 10 * LDEBetaValue
            startCell.offset(currentRow - startCell.Row + 10 - i - 1, LDEColumn.Column - startCell.Column).Value = mixture
        Next i

        ' Move the output row down for the next photon value's calculations
        currentRow = currentRow + 10
    Next photonCell

    ' Reset the currentRow for SDE calculations
    currentRow = startCell.Row
    ' SDE Mixture Calculations
    For Each photonCell In photonCells
        photon = photonCell.Value
        ' Use the row of photonCell to fetch the SDE photon value
        SDEPhotonValue = SDEColumn.Cells(photonCell.Row - photonCells.Cells(1).Row + 1, 1).Value
        ' Fetch the SDE beta value
        SDEBetaValue = betaCell.offset(0, SDEColumn.Column - startCell.Column).Value

        ' Calculate the mixtures for the SDE column
        For i = 9 To 1 Step -1
            mixture = (10 - i) / 10 * SDEPhotonValue + i / 10 * SDEBetaValue
            startCell.offset(currentRow - startCell.Row + 10 - i - 1, SDEColumn.Column - startCell.Column).Value = mixture
        Next i

        ' Move the output row down for the next photon value's calculations
        currentRow = currentRow + 10
    Next photonCell

    ' Reset the currentRow for Beta SDE calculations
    currentRow = startCell.Row
    ' Beta SDE Mixture Calculations
    For Each photonCell In photonCells
        photon = photonCell.Value
        ' Use the row of photonCell to fetch the Beta SDE photon value
        BetaSDEPhotonValue = BetaSDEColumn.Cells(photonCell.Row - photonCells.Cells(1).Row + 1, 1).Value
        ' Fetch the Beta SDE beta value
        BetaSDEBetaValue = betaCell.offset(0, BetaSDEColumn.Column - startCell.Column).Value

        ' Calculate the mixtures for the Beta SDE column
        For i = 9 To 1 Step -1
            mixture = (10 - i) / 10 * BetaSDEPhotonValue + i / 10 * BetaSDEBetaValue
            startCell.offset(currentRow - startCell.Row + 10 - i - 1, BetaSDEColumn.Column - startCell.Column).Value = mixture
        Next i

        ' Move the output row down for the next photon value's calculations
        currentRow = currentRow + 10
    Next photonCell
End Sub

'Sub PopulateMixturesForBetaAndPhotons()
'    Dim startCell As Range
'    Dim betaCell As Range
'    Dim betaNameCell As Range
'    Dim photonCells As Range
'    Dim photonCell As Range
'    Dim photonNameColumn As Range
'    Dim DDEColumn As Range
'    Dim LDEColumn As Range
'    Dim SDEColumn As Range
'    Dim BetaSDEColumn As Range
'    Dim beta As Double
'    Dim betaName As String
'    Dim i As Integer, k As Integer
'    Dim mixture As Double
'    Dim ratio As String
'    Dim currentRow As Long
'    Dim photon As Double
'    Dim photonType As String
'    Dim DDEPhotonValue As Double
'    Dim LDEPhotonValue As Double
'    Dim SDEPhotonValue As Double
'    Dim BetaSDEPhotonValue As Double
'    Dim DDEBetaValue As Double
'    Dim LDEBetaValue As Double
'    Dim SDEBetaValue As Double
'    Dim BetaSDEBetaValue As Double
'
'    ' Prompt user to select the starting cell for numeric output
'    On Error Resume Next
'    Set startCell = Application.InputBox("Select the starting cell for numeric output:", Type:=8)
'    On Error GoTo 0
'    If startCell Is Nothing Then Exit Sub ' User canceled
'
'    ' Prompt user to select the cell containing the beta value
'    On Error Resume Next
'    Set betaCell = Application.InputBox("Select the cell containing the beta value:", Type:=8)
'    On Error GoTo 0
'    If betaCell Is Nothing Then Exit Sub ' User canceled
'
'    ' Prompt user to select the cell containing the beta name
'    On Error Resume Next
'    Set betaNameCell = Application.InputBox("Select the cell containing the beta name:", Type:=8)
'    On Error GoTo 0
'    If betaNameCell Is Nothing Then Exit Sub ' User canceled
'
'    ' Get the beta value and name
'    beta = betaCell.Value
'    betaName = betaNameCell.Value
'
'    ' Prompt user to select the range of photon values (allowing non-contiguous selection)
'    On Error Resume Next
'    Set photonCells = Application.InputBox("Select the cells containing the photon values (you can select multiple non-contiguous cells):", Type:=8)
'    On Error GoTo 0
'    If photonCells Is Nothing Then Exit Sub ' User canceled
'
'    ' Prompt user to select the column containing the photon names
'    On Error Resume Next
'    Set photonNameColumn = Application.InputBox("Select the column containing the photon names:", Type:=8)
'    On Error GoTo 0
'    If photonNameColumn Is Nothing Then Exit Sub ' User canceled
'
'    ' Prompt user to select the column for DDE data
'    On Error Resume Next
'    Set DDEColumn = Application.InputBox("Select the column containing the DDE data:", Type:=8)
'    On Error GoTo 0
'    If DDEColumn Is Nothing Then Exit Sub ' User canceled
'
'    ' Prompt user to select the column for LDE data
'    On Error Resume Next
'    Set LDEColumn = Application.InputBox("Select the column containing the LDE data:", Type:=8)
'    On Error GoTo 0
'    If LDEColumn Is Nothing Then Exit Sub ' User canceled
'
'    ' Prompt user to select the column for SDE data
'    On Error Resume Next
'    Set SDEColumn = Application.InputBox("Select the column containing the SDE data:", Type:=8)
'    On Error GoTo 0
'    If SDEColumn Is Nothing Then Exit Sub ' User canceled
'
'    ' Prompt user to select the column for Beta SDE data
'    On Error Resume Next
'    Set BetaSDEColumn = Application.InputBox("Select the column containing the Beta SDE data:", Type:=8)
'    On Error GoTo 0
'    If BetaSDEColumn Is Nothing Then Exit Sub ' User canceled
'
'    ' Determine the start row for calculations
'    currentRow = startCell.Row
'
'    ' First, print the ratios and initial mixtures for the first set of columns
'    For Each photonCell In photonCells
'        photon = photonCell.Value
'        ' Use the row of photonCell to fetch the photon name
'        photonType = photonNameColumn.Cells(photonCell.Row - photonCells.Cells(1).Row + 1, 1).Value
'
'        ' Populate the columns with types and ratios for this set of photon values (once)
'        For i = 1 To 9
'            ratio = i & ":" & (10 - i)
'            startCell.Offset(currentRow - startCell.Row + i - 1, -3).Value = photonType
'            startCell.Offset(currentRow - startCell.Row + i - 1, -2).Value = betaName
'            startCell.Offset(currentRow - startCell.Row + i - 1, -1).Value = ratio
'        Next i
'
'        ' Calculate the mixtures for the first photon value
'        For i = 9 To 1 Step -1
'            mixture = (10 - i) / 10 * photon + i / 10 * beta
'            startCell.Offset(currentRow - startCell.Row + 10 - i - 1, 0).Value = mixture
'        Next i
'
'        ' Move the output row down for the next photon value's calculations
'        currentRow = currentRow + 10 ' 9 rows of mixtures + 1 row for spacing
'    Next photonCell
'
'    ' Reset the currentRow for additional columns of mixtures
'    currentRow = startCell.Row
'
'    ' Loop through the next 3 sets of photon and beta values (columns)
'    For k = 1 To 3 ' 1 is the second column, 2 is the third, etc.
'        ' Adjust for the current column offset (starting from the second set of columns)
'        colOffset = k
'        calcColOffset = colOffset * 1 ' Adjust the number of columns to move to the right for each set
'
'        ' Get the current beta and beta name from the correct column
'        beta = betaCell.Offset(0, colOffset).Value
'        betaName = betaNameCell.Offset(0, colOffset).Value
'
'        ' Loop through each selected photon value and corresponding type
'        For Each photonCell In photonCells
'            photon = photonCell.Offset(0, colOffset).Value
'
'            ' Calculate the mixtures for the current photon value
'            For i = 9 To 1 Step -1
'                mixture = (10 - i) / 10 * photon + i / 10 * beta
'                startCell.Offset(currentRow - startCell.Row + 10 - i - 1, calcColOffset).Value = mixture
'            Next i
'
'            ' Move the output row down for the next photon value's calculations
'            currentRow = currentRow + 10 ' 9 rows of mixtures + 1 row for spacing
'        Next photonCell
'
'        ' Reset the currentRow for the next set of columns
'        currentRow = startCell.Row
'    Next k
'
'    ' DDE Mixture Calculations
'    currentRow = startCell.Row
'    For Each photonCell In photonCells
'        photon = photonCell.Value
'        ' Use the row of photonCell to fetch the DDE photon value
'        DDEPhotonValue = DDEColumn.Cells(photonCell.Row - photonCells.Cells(1).Row + 1, 1).Value
'        ' Fetch the DDE beta value
'        DDEBetaValue = betaCell.Offset(0, DDEColumn.Column - startCell.Column).Value
'
'        ' Calculate the mixtures for the DDE column
'        For i = 9 To 1 Step -1
'            mixture = (10 - i) / 10 * DDEPhotonValue + i / 10 * DDEBetaValue
'            startCell.Offset(currentRow - startCell.Row + 10 - i - 1, DDEColumn.Column - startCell.Column).Value = mixture
'        Next i
'
'        ' Move the output row down for the next photon value's calculations
'        currentRow = currentRow + 10
'    Next photonCell
'
'    ' Reset the currentRow for LDE calculations
'    currentRow = startCell.Row
'    ' LDE Mixture Calculations
'    For Each photonCell In photonCells
'        photon = photonCell.Value
'        ' Use the row of photonCell to fetch the LDE photon value
'        LDEPhotonValue = LDEColumn.Cells(photonCell.Row - photonCells.Cells(1).Row + 1, 1).Value
'        ' Fetch the LDE beta value
'        LDEBetaValue = betaCell.Offset(0, LDEColumn.Column - startCell.Column).Value
'
'        ' Calculate the mixtures for the LDE column
'        For i = 9 To 1 Step -1
'            mixture = (10 - i) / 10 * LDEPhotonValue + i / 10 * LDEBetaValue
'            startCell.Offset(currentRow - startCell.Row + 10 - i - 1, LDEColumn.Column - startCell.Column).Value = mixture
'        Next i
'
'        ' Move the output row down for the next photon value's calculations
'        currentRow = currentRow + 10
'    Next photonCell
'
'    ' Reset the currentRow for SDE calculations
'    currentRow = startCell.Row
'    ' SDE Mixture Calculations
'    For Each photonCell In photonCells
'        photon = photonCell.Value
'        ' Use the row of photonCell to fetch the SDE photon value
'        SDEPhotonValue = SDEColumn.Cells(photonCell.Row - photonCells.Cells(1).Row + 1, 1).Value
'        ' Fetch the SDE beta value
'        SDEBetaValue = betaCell.Offset(0, SDEColumn.Column - startCell.Column).Value
'
'        ' Calculate the mixtures for the SDE column
'        For i = 9 To 1 Step -1
'            mixture = (10 - i) / 10 * SDEPhotonValue + i / 10 * SDEBetaValue
'            startCell.Offset(currentRow - startCell.Row + 10 - i - 1, SDEColumn.Column - startCell.Column).Value = mixture
'        Next i
'
'        ' Move the output row down for the next photon value's calculations
'        currentRow = currentRow + 10
'    Next photonCell
'
'    ' Reset the currentRow for Beta SDE calculations
'    currentRow = startCell.Row
'    ' Beta SDE Mixture Calculations
'    For Each photonCell In photonCells
'        photon = photonCell.Value
'        ' Use the row of photonCell to fetch the Beta SDE photon value
'        BetaSDEPhotonValue = BetaSDEColumn.Cells(photonCell.Row - photonCells.Cells(1).Row + 1, 1).Value
'        ' Fetch the Beta SDE beta value
'        BetaSDEBetaValue = betaCell.Offset(0, BetaSDEColumn.Column - startCell.Column).Value
'
'        ' Calculate the mixtures for the Beta SDE column
'        For i = 9 To 1 Step -1
'            mixture = (10 - i) / 10 * BetaSDEPhotonValue + i / 10 * BetaSDEBetaValue
'            startCell.Offset(currentRow - startCell.Row + 10 - i - 1, BetaSDEColumn.Column - startCell.Column).Value = mixture
'        Next i
'
'        ' Move the output row down for the next photon value's calculations
'        currentRow = currentRow + 10
'    Next photonCell
'End Sub

Sub PopulateOldAlgorithmMixturesForBetaAndPhotons()
    Dim startCell As Range
    Dim betaCell As Range
    Dim betaNameCell As Range
    Dim photonCells As Range
    Dim photonCell As Range
    Dim photonNameColumn As Range
    Dim LDEColumn As Range
    Dim SDEColumn As Range
    Dim BetaSDEColumn As Range
    Dim AlgDDEColumn As Range
    Dim AlgSDEColumn As Range
    Dim beta As Double
    Dim betaName As String
    Dim i As Integer, k As Integer
    Dim mixture As Double
    Dim ratio As String
    Dim currentRow As Long
    Dim photon As Double
    Dim photonType As String
    Dim LDEPhotonValue As Double
    Dim SDEPhotonValue As Double
    Dim BetaSDEPhotonValue As Double
    Dim AlgDDEPhotonValue As Double
    Dim AlgSDEPhotonValue As Double
    Dim LDEBetaValue As Double
    Dim SDEBetaValue As Double
    Dim BetaSDEBetaValue As Double
    Dim AlgDDEBetaValue As Double
    Dim AlgSDEBetaValue As Double

    ' Prompt user to select the starting cell for numeric output
    On Error Resume Next
    Set startCell = Application.InputBox("Select the starting cell for numeric output:", Type:=8)
    On Error GoTo 0
    If startCell Is Nothing Then Exit Sub ' User canceled

    ' Prompt user to select the cell containing the beta value
    On Error Resume Next
    Set betaCell = Application.InputBox("Select the cell containing the beta value:", Type:=8)
    On Error GoTo 0
    If betaCell Is Nothing Then Exit Sub ' User canceled

    ' Prompt user to select the cell containing the beta name
    On Error Resume Next
    Set betaNameCell = Application.InputBox("Select the cell containing the beta name:", Type:=8)
    On Error GoTo 0
    If betaNameCell Is Nothing Then Exit Sub ' User canceled

    ' Get the beta value and name
    beta = betaCell.Value
    betaName = betaNameCell.Value

    ' Prompt user to select the range of photon values (allowing non-contiguous selection)
    On Error Resume Next
    Set photonCells = Application.InputBox("Select the cells containing the photon values (you can select multiple non-contiguous cells):", Type:=8)
    On Error GoTo 0
    If photonCells Is Nothing Then Exit Sub ' User canceled

    ' Prompt user to select the column containing the photon names
    On Error Resume Next
    Set photonNameColumn = Application.InputBox("Select the column containing the photon names:", Type:=8)
    On Error GoTo 0
    If photonNameColumn Is Nothing Then Exit Sub ' User canceled

    ' Prompt user to select the column for LDE data
    On Error Resume Next
    Set LDEColumn = Application.InputBox("Select the column containing the LDE data:", Type:=8)
    On Error GoTo 0
    If LDEColumn Is Nothing Then Exit Sub ' User canceled

    ' Prompt user to select the column for SDE data
    On Error Resume Next
    Set SDEColumn = Application.InputBox("Select the column containing the SDE data:", Type:=8)
    On Error GoTo 0
    If SDEColumn Is Nothing Then Exit Sub ' User canceled

    ' Prompt user to select the column for Beta SDE data
    On Error Resume Next
    Set BetaSDEColumn = Application.InputBox("Select the column containing the Beta SDE data:", Type:=8)
    On Error GoTo 0
    If BetaSDEColumn Is Nothing Then Exit Sub ' User canceled

    ' Prompt user to select the column for AlgDDE data
    On Error Resume Next
    Set AlgDDEColumn = Application.InputBox("Select the column containing the AlgDDE data:", Type:=8)
    On Error GoTo 0
    If AlgDDEColumn Is Nothing Then Exit Sub ' User canceled

    ' Prompt user to select the column for AlgSDE data
    On Error Resume Next
    Set AlgSDEColumn = Application.InputBox("Select the column containing the AlgSDE data:", Type:=8)
    On Error GoTo 0
    If AlgSDEColumn Is Nothing Then Exit Sub ' User canceled

    ' Determine the start row for calculations
    currentRow = startCell.Row

    ' First, print the ratios and initial mixtures for the first set of columns
    For Each photonCell In photonCells
        photon = photonCell.Value
        ' Use the row of photonCell to fetch the photon name
        photonType = photonNameColumn.Cells(photonCell.Row - photonCells.Cells(1).Row + 1, 1).Value

        ' Populate the columns with types and ratios for this set of photon values (once)
        For i = 1 To 9
            ratio = i & ":" & (10 - i)
            startCell.offset(currentRow - startCell.Row + i - 1, -3).Value = photonType
            startCell.offset(currentRow - startCell.Row + i - 1, -2).Value = betaName
            startCell.offset(currentRow - startCell.Row + i - 1, -1).Value = ratio
        Next i

        ' Calculate the mixtures for the first photon value
        For i = 9 To 1 Step -1
            mixture = (10 - i) / 10 * photon + i / 10 * beta
            startCell.offset(currentRow - startCell.Row + 10 - i - 1, 0).Value = mixture
        Next i

        ' Move the output row down for the next photon value's calculations
        currentRow = currentRow + 10 ' 9 rows of mixtures + 1 row for spacing
    Next photonCell

    ' Reset the currentRow for additional columns of mixtures
    currentRow = startCell.Row

    ' Loop through the next 4 sets of photon and beta values (columns)
    For k = 1 To 4 ' 1 is the second column, 2 is the third, etc.
        ' Adjust for the current column offset (starting from the second set of columns)
        colOffset = k
        calcColOffset = colOffset * 1 ' Adjust the number of columns to move to the right for each set

        ' Get the current beta and beta name from the correct column
        beta = betaCell.offset(0, colOffset).Value
        betaName = betaNameCell.offset(0, colOffset).Value

        ' Loop through each selected photon value and corresponding type
        For Each photonCell In photonCells
            photon = photonCell.offset(0, colOffset).Value

            ' Calculate the mixtures for the current photon value
            For i = 9 To 1 Step -1
                mixture = (10 - i) / 10 * photon + i / 10 * beta
                startCell.offset(currentRow - startCell.Row + 10 - i - 1, calcColOffset).Value = mixture
            Next i

            ' Move the output row down for the next photon value's calculations
            currentRow = currentRow + 10 ' 9 rows of mixtures + 1 row for spacing
        Next photonCell

        ' Reset the currentRow for the next set of columns
        currentRow = startCell.Row
    Next k

    ' Calculate LDE mixtures
    currentRow = startCell.Row
    For Each photonCell In photonCells
        photon = photonCell.Value
        ' Use the row of photonCell to fetch the LDE photon value
        LDEPhotonValue = LDEColumn.Cells(photonCell.Row - photonCells.Cells(1).Row + 1, 1).Value
        ' Fetch the LDE beta value
        LDEBetaValue = betaCell.offset(0, LDEColumn.Column - startCell.Column).Value

        ' Calculate the mixtures for the LDE column
        For i = 9 To 1 Step -1
            mixture = (10 - i) / 10 * LDEPhotonValue + i / 10 * LDEBetaValue
            startCell.offset(currentRow - startCell.Row + 10 - i - 1, LDEColumn.Column - startCell.Column).Value = mixture
        Next i

        ' Move the output row down for the next photon value's calculations
        currentRow = currentRow + 10
    Next photonCell

    ' Reset the currentRow for SDE calculations
    currentRow = startCell.Row
    ' SDE Mixture Calculations
    For Each photonCell In photonCells
        photon = photonCell.Value
        ' Use the row of photonCell to fetch the SDE photon value
        SDEPhotonValue = SDEColumn.Cells(photonCell.Row - photonCells.Cells(1).Row + 1, 1).Value
        ' Fetch the SDE beta value
        SDEBetaValue = betaCell.offset(0, SDEColumn.Column - startCell.Column).Value

        ' Calculate the mixtures for the SDE column
        For i = 9 To 1 Step -1
            mixture = (10 - i) / 10 * SDEPhotonValue + i / 10 * SDEBetaValue
            startCell.offset(currentRow - startCell.Row + 10 - i - 1, SDEColumn.Column - startCell.Column).Value = mixture
        Next i

        ' Move the output row down for the next photon value's calculations
        currentRow = currentRow + 10
    Next photonCell

    ' Reset the currentRow for Beta SDE calculations
    currentRow = startCell.Row
    ' Beta SDE Mixture Calculations
    For Each photonCell In photonCells
        photon = photonCell.Value
        ' Use the row of photonCell to fetch the Beta SDE photon value
        BetaSDEPhotonValue = BetaSDEColumn.Cells(photonCell.Row - photonCells.Cells(1).Row + 1, 1).Value
        ' Fetch the Beta SDE beta value
        BetaSDEBetaValue = betaCell.offset(0, BetaSDEColumn.Column - startCell.Column).Value

        ' Calculate the mixtures for the Beta SDE column
        For i = 9 To 1 Step -1
            mixture = (10 - i) / 10 * BetaSDEPhotonValue + i / 10 * BetaSDEBetaValue
            startCell.offset(currentRow - startCell.Row + 10 - i - 1, BetaSDEColumn.Column - startCell.Column).Value = mixture
        Next i

        ' Move the output row down for the next photon value's calculations
        currentRow = currentRow + 10
    Next photonCell

    ' Reset the currentRow for AlgDDE calculations
    currentRow = startCell.Row
    ' AlgDDE Mixture Calculations
    For Each photonCell In photonCells
        photon = photonCell.Value
        ' Use the row of photonCell to fetch the AlgDDE photon value
        AlgDDEPhotonValue = AlgDDEColumn.Cells(photonCell.Row - photonCells.Cells(1).Row + 1, 1).Value
        ' Fetch the AlgDDE beta value
        AlgDDEBetaValue = betaCell.offset(0, AlgDDEColumn.Column - startCell.Column).Value

        ' Calculate the mixtures for the AlgDDE column
        For i = 9 To 1 Step -1
            mixture = (10 - i) / 10 * AlgDDEPhotonValue + i / 10 * AlgDDEBetaValue
            startCell.offset(currentRow - startCell.Row + 10 - i - 1, AlgDDEColumn.Column - startCell.Column).Value = mixture
        Next i

        ' Move the output row down for the next photon value's calculations
        currentRow = currentRow + 10
    Next photonCell

    ' Reset the currentRow for AlgSDE calculations
    currentRow = startCell.Row
    ' AlgSDE Mixture Calculations
    For Each photonCell In photonCells
        photon = photonCell.Value
        ' Use the row of photonCell to fetch the AlgSDE photon value
        AlgSDEPhotonValue = AlgSDEColumn.Cells(photonCell.Row - photonCells.Cells(1).Row + 1, 1).Value
        ' Fetch the AlgSDE beta value
        AlgSDEBetaValue = betaCell.offset(0, AlgSDEColumn.Column - startCell.Column).Value

        ' Calculate the mixtures for the AlgSDE column
        For i = 9 To 1 Step -1
            mixture = (10 - i) / 10 * AlgSDEPhotonValue + i / 10 * AlgSDEBetaValue
            startCell.offset(currentRow - startCell.Row + 10 - i - 1, AlgSDEColumn.Column - startCell.Column).Value = mixture
        Next i

        ' Move the output row down for the next photon value's calculations
        currentRow = currentRow + 10
    Next photonCell

End Sub


Function CalculateMeanPercentDifference(rng1 As Range, rng2 As Range) As Double
    Dim i As Long
    Dim sumPercentDiff As Double
    Dim count As Long
    Dim percentDiff As Double

    ' Initialize variables
    sumPercentDiff = 0
    count = 0

    ' Ensure both ranges have the same number of cells
    If rng1.Cells.count <> rng2.Cells.count Then
        CalculateMeanPercentDifference = CVErr(xlErrNA)
        Exit Function
    End If

    ' Loop through each cell in the ranges
    For i = 1 To rng1.Cells.count
        ' Check if both cells are numeric and not empty or an error
        If IsNumeric(rng1.Cells(i).Value) And IsNumeric(rng2.Cells(i).Value) And _
           Not IsEmpty(rng1.Cells(i).Value) And Not IsEmpty(rng2.Cells(i).Value) And _
           Not IsError(rng1.Cells(i).Value) And Not IsError(rng2.Cells(i).Value) Then

            If rng1.Cells(i).Value <> 0 Then
                percentDiff = Abs((rng1.Cells(i).Value - rng2.Cells(i).Value) / rng1.Cells(i).Value) * 100
                sumPercentDiff = sumPercentDiff + percentDiff
                count = count + 1
            End If
        End If
    Next i

    ' Calculate Mean Percent Difference
    If count > 0 Then
        CalculateMeanPercentDifference = sumPercentDiff / count
    Else
        CalculateMeanPercentDifference = CVErr(xlErrDiv0)
    End If
End Function

Function CalculateRMSE(rng1 As Range, rng2 As Range) As Double
    Dim i As Long
    Dim sumSqErrors As Double
    Dim count As Long
    Dim diff As Double

    ' Initialize variables
    sumSqErrors = 0
    count = 0

    ' Ensure both ranges have the same number of cells
    If rng1.Cells.count <> rng2.Cells.count Then
        CalculateRMSE = CVErr(xlErrNA)
        Exit Function
    End If

    ' Loop through each cell in the ranges
    For i = 1 To rng1.Cells.count
        ' Check if both cells are numeric and not empty or an error
        If IsNumeric(rng1.Cells(i).Value) And IsNumeric(rng2.Cells(i).Value) And _
           Not IsEmpty(rng1.Cells(i).Value) And Not IsEmpty(rng2.Cells(i).Value) And _
           Not IsError(rng1.Cells(i).Value) And Not IsError(rng2.Cells(i).Value) Then

            diff = rng1.Cells(i).Value - rng2.Cells(i).Value
            sumSqErrors = sumSqErrors + diff ^ 2
            count = count + 1
        End If
    Next i

    ' Calculate RMSE
    If count > 0 Then
        CalculateRMSE = Sqr(sumSqErrors / count)
    Else
        CalculateRMSE = CVErr(xlErrDiv0)
    End If
End Function

Public Function CalculatePercentDifference(experimental_value As Double, known_value As Double) As Variant
    Dim percentDifference As Double
    
    ' Check if both values are valid numeric entries and not errors or empty
    If Not IsEmpty(known_value) And Not IsEmpty(experimental_value) And _
       Not IsError(known_value) And IsNumeric(known_value) And _
       Not IsError(experimental_value) And IsNumeric(experimental_value) Then

        ' Calculate percent difference
        If known_value <> 0 Then
            percentDifference = Abs(experimental_value - known_value) / Abs(known_value) * 100
            CalculatePercentDifference = percentDifference
        Else
            CalculatePercentDifference = CVErr(xlErrDiv0)
        End If
    Else
        CalculatePercentDifference = ""
    End If
End Function


Public Function CalculateSignedPercentDifference(experimental_value As Double, known_value As Double) As Variant
    Dim percentDifference As Double
    
    ' Check if both values are valid numeric entries and not errors or empty
    If Not IsEmpty(known_value) And Not IsEmpty(experimental_value) And _
       Not IsError(known_value) And IsNumeric(known_value) And _
       Not IsError(experimental_value) And IsNumeric(experimental_value) Then

        ' Calculate percent difference
        If known_value <> 0 Then
            percentDifference = (experimental_value - known_value) / Abs(known_value) * 100
            CalculateSignedPercentDifference = percentDifference
        Else
            CalculateSignedPercentDifference = CVErr(xlErrDiv0)
        End If
    Else
        CalculateSignedPercentDifference = ""
    End If
End Function

Function determine_within_asym_tolerances(known_value As Variant, experimental_value As Variant, lower_tolerance As Double, upper_tolerance As Double) As Variant
    Dim lower_limit As Double
    Dim upper_limit As Double
    Dim deviation As Double
    
    ' Check if known_value is empty, zero, or not a number
    If IsEmpty(known_value) Or IsError(known_value) Or Not IsNumeric(known_value) Or known_value = 0 Then
        determine_within_asym_tolerances = "-"  ' Return blank if known_value is invalid
        Exit Function
    End If
    
    ' Calculate the tolerance limits
    lower_limit = known_value * ((100 - lower_tolerance) / 100)
    upper_limit = known_value * ((100 + upper_tolerance) / 100)
    
    ' Check if experimental_value is within the tolerance range
    If experimental_value >= lower_limit And experimental_value <= upper_limit Then
        determine_within_asym_tolerances = 1  ' Return 1 if within range
    Else
        ' Calculate deviation in percentage
        If experimental_value < lower_limit Then
            deviation = -((known_value - experimental_value) / known_value) * 100
        Else
            deviation = ((experimental_value - known_value) / known_value) * 100
        End If
        
        determine_within_asym_tolerances = Round(deviation, 2)  ' Return the deviation, rounded to 2 decimal places
    End If
End Function


Function determine_within_tolerances(known_value As Variant, experimental_value As Variant, tolerance As Double) As Variant
    Dim lower_limit As Double
    Dim upper_limit As Double
    Dim deviation As Double
    
    ' Check if known_value is empty, zero, or not a number
    If IsEmpty(known_value) Or IsError(known_value) Or Not IsNumeric(known_value) Or known_value = 0 Then
        determine_within_tolerances = ""  ' Return blank if known_value is invalid
        Exit Function
    End If
    
    ' Calculate the tolerance limits
    lower_limit = known_value * ((100 - tolerance) / 100)
    upper_limit = known_value * ((100 + tolerance) / 100)
    
    ' Check if experimental_value is within the tolerance range
    If experimental_value >= lower_limit And experimental_value <= upper_limit Then
        determine_within_tolerances = 1  ' Return 1 if within range
    Else
        ' Calculate deviation in percentage
        If experimental_value < lower_limit Then
            deviation = -((known_value - experimental_value) / known_value) * 100
        Else
            deviation = ((experimental_value - known_value) / known_value) * 100
        End If
        
        determine_within_tolerances = Round(deviation, 2)  ' Return the deviation, rounded to 2 decimal places
    End If
End Function

Function calc_beta_pass(experimental_beta As Double, known_beta As Double)
    Dim beta_pass As Integer
    
    If (experimental_beta = known_beta) Then
        beta_pass = 1
    ElseIf (experimental_beta = 1) And (known_beta > 2) Then
        beta_pass = 1
    Else
        beta_pass = 0
    End If
    
    calc_beta_pass = beta_pass

End Function

Function calc_pass_fail_pct(inputs As Range, var_type As String) As Double
    Dim NumberOfCells As Double
    Dim NumberOfPasses As Double
    Dim NumberOfValidCells As Double
    
    NumberOfCells = inputs.Cells.count
    NumberOfPasses = 0
    NumberOfValidCells = 0
    
    Dim pct_wrong As Double
    
    For Each cell In inputs
        If Not IsEmpty(cell.Value) And Not IsError(cell.Value) And IsNumeric(cell.Value) Then
            NumberOfValidCells = NumberOfValidCells + 1
            If cell.Value = var_type Then
                NumberOfPasses = NumberOfPasses + 1
            End If
        End If
    Next cell
    
    If NumberOfValidCells > 0 Then
        pct_wrong = (NumberOfValidCells - NumberOfPasses) / NumberOfValidCells
    Else
        pct_wrong = 0 ' Handle the case where there are no valid cells
    End If
    
    calc_pass_fail_pct = pct_wrong
    
End Function

' Function to generate a Gaussian value based on the given amplitude, mean, and standard deviation
Function GaussianValue(A As Double, mu As Double, sigma As Double, x As Double) As Double
    GaussianValue = A * Exp(-0.5 * ((x - mu) / sigma) ^ 2)
End Function

' Subroutine to add Gaussian noise to the first four cells in each row of the specified range
Sub AddGaussianNoiseToFirstFourCellsByRow(cellRange As Range, amplitudePercent As Double, printColumn As Range)
    Dim rowRange As Range
    Dim noise As Double
    Dim amplitude As Double
    Dim lowestCellValue As Double
    Dim i As Integer
    Dim SkipRow As Boolean
    Dim rowNum As Long
    Dim mu As Double
    Dim sigma As Double
    Dim x As Double
    Dim fnX As Double
    
    ' Loop through each row in the specified range
    For Each rowRange In cellRange.Rows
        ' Initialize skipRow flag
        SkipRow = False
        
        ' Check if any of the first four cells in the row are empty
        For i = 1 To 4
            If IsEmpty(rowRange.Cells(1, i).Value) Then
                SkipRow = True
                Exit For
            End If
        Next i
        
        ' If the row contains an empty cell, skip to the next row
        If SkipRow Then GoTo SkipRow
        
        ' Find the lowest value among the first four cells in the row (A)
        lowestCellValue = rowRange.Cells(1, 1).Value
        For i = 2 To 4
            If rowRange.Cells(1, i).Value < lowestCellValue Then
                lowestCellValue = rowRange.Cells(1, i).Value
            End If
        Next i
        
        ' Define Gaussian parameters
        amplitude = lowestCellValue ' A is both the amplitude and the mean (mu)
        mu = amplitude
        sigma = 0
        sigma = amplitude * amplitudePercent / 100 ' sigma = A * (user-defined percentage)
        
        ' Generate a random x within the range (mu - 4*sigma) < x < (mu + 4*sigma)
        x = 0
        x = mu - 4 * sigma + 8 * sigma * Rnd
        
        ' Evaluate the Gaussian at the random x
        fnX = GaussianValue(amplitude, mu, sigma, x)
        
        ' Calculate the noise as the absolute difference between A and the Gaussian value at x
        noise = Abs(amplitude - x)
        
        ' Apply the same noise value to the first four cells in the row
        For i = 1 To 4
            rowRange.Cells(1, i).Value = rowRange.Cells(1, i).Value + noise
        Next i
        
        ' Print the noise value in the specified column of the same row
        printColumn.Cells(rowRange.Row - 1, 1).Value = noise
SkipRow:
    Next rowRange
End Sub

' Main subroutine to prompt for cell range, amplitude percentage, and print column
Sub AddGaussianNoise()
    Dim cellRange As Range
    Dim amplitudeInput As Variant
    Dim amplitudePercent As Double
    Dim prompt As String
    Dim printColumn As Range
    
    ' Prompt the user to select a range of cells
    On Error Resume Next
    Set cellRange = Application.InputBox("Select the range of cells to apply Gaussian noise:", Type:=8)
    On Error GoTo 0
    
    ' If no range was selected, exit the subroutine
    If cellRange Is Nothing Then
        MsgBox "No range selected. Operation canceled."
        Exit Sub
    End If
    
    ' Prompt the user to select a column where the noise values will be printed
    On Error Resume Next
    Set printColumn = Application.InputBox("Select the column to print the noise values:", Type:=8)
    On Error GoTo 0
    
    ' If no column was selected, exit the subroutine
    If printColumn Is Nothing Then
        MsgBox "No column selected. Operation canceled."
        Exit Sub
    End If
    
    ' Prompt the user to enter the amplitude percentage
    prompt = "Enter the amplitude percentage (e.g., 10 for 10%):"
    amplitudeInput = InputBox(prompt, "Set Amplitude Percentage", 10)
    
    ' Validate the amplitude percentage input
    If Not IsNumeric(amplitudeInput) Or Val(amplitudeInput) <= 0 Then
        MsgBox "Invalid amplitude percentage entered. Operation canceled."
        Exit Sub
    End If
    
    ' Explicitly convert amplitudeInput to a Double
    amplitudePercent = CDbl(amplitudeInput)
    
    ' Call the function to add Gaussian noise to the selected cells by row
    AddGaussianNoiseToFirstFourCellsByRow cellRange, amplitudePercent, printColumn
End Sub

