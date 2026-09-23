$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$exportDir = Join-Path $projectRoot 'reports/presentation'
$deckPath = Join-Path $exportDir 'GWAS_First_Lab_Presentation_CRC_completed.pptx'
$existingPids = @(Get-Process POWERPNT -ErrorAction SilentlyContinue | ForEach-Object Id)
$app = $null
$deck = $null
try {
    $app = New-Object -ComObject PowerPoint.Application
    $deck = $app.Presentations.Open($deckPath, -1, 0, 0)
    if ($deck.Slides.Count -ne 10) { throw 'Expected ten slides' }
    $deck.Slides.Item(9).Export((Join-Path $exportDir 'slide09_reference.png'), 'PNG', 1920, 1080)
    $deck.Slides.Item(10).Export((Join-Path $exportDir 'slide10_preview.png'), 'PNG', 1920, 1080)
    $check = @()
    foreach ($s in $deck.Slides.Item(10).Shapes) {
        if ($s.HasTextFrame -and $s.TextFrame2.HasText) {
            $t = $s.TextFrame2.TextRange
            $check += [pscustomobject]@{
                Name=$s.Name; Text=$t.Text; Width=$s.Width; Height=$s.Height
                TextWidth=$t.BoundWidth; TextHeight=$t.BoundHeight
                Overflow=($t.BoundWidth -gt ($s.Width + 1) -or $t.BoundHeight -gt ($s.Height + 1))
            }
        }
    }
    $check | ConvertTo-Json -Depth 3 | Set-Content -LiteralPath (Join-Path $exportDir 'powerpoint_layout_check.json') -Encoding UTF8
    $overflow = @($check | Where-Object Overflow)
    Write-Output "PowerPoint opened ten-slide deck. Text overflow boxes: $($overflow.Count)"
    $overflow | Format-List
    if ($overflow.Count -gt 0) { throw 'Text overflow detected; adjust slide layout' }
} finally {
    if ($null -ne $deck) { $deck.Close(); [void][System.Runtime.InteropServices.Marshal]::ReleaseComObject($deck) }
    if ($null -ne $app) {
        $newPids = @(Get-Process POWERPNT -ErrorAction SilentlyContinue | Where-Object { $_.Id -notin $existingPids })
        if ($existingPids.Count -eq 0 -and $newPids.Count -gt 0 -and $app.Presentations.Count -eq 0) { $app.Quit() }
        [void][System.Runtime.InteropServices.Marshal]::ReleaseComObject($app)
    }
}
