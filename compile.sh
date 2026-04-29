cwd=$(pwd)
cd "$cwd/playground/CMSSW_15_0_2/src/HeavyFlavorAnalysor/Onia2MuMu"
cp "$cwd"/src/* src/
cp "$cwd"/interface/* interface/
scramv1 b && echo "succeeded" || echo "failed"
cd "$cwd"
