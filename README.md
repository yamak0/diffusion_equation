## 三相脳循環モデル

# このコードは以下のライブラリに依存しています．
1. cmake (https://cmake.org/)
2. HDF5 (https://www.hdfgroup.org/solutions/hdf5/)
3. intel compiler (oneAPI) (https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#gs.enmkd3)
4. TextParser(https://github.com/avr-aics-riken/TextParser)

# 必要なinput file
- node.dat : 計算領域の節点座標
- element.dat : 有限要素の構成情報
- boundary_fluid.dat : 脳脊髄液領域の境界条件
- fluid_phi.dat : 各要素の脳脊髄液領域の相割合
- fluid_node_phi.dat : 各要素の脳脊髄液領域の相割合を節点に変換したもの
- boundary_solid.dat : 固体領域の境界条件
- solid_phi.dat : 各要素の固体領域の相割合
- solid_node_phi.dat : 各要素の固体領域の相割合を節点に変換したもの
- boundary_vessel.dat : 血管領域の境界条件
- vessel_phi.dat : 各要素の血管領域の相割合
- vessel_node_phi.dat : 各要素の血管領域の相割合を節点に変換したもの

# コードの動かし方
build.shを実行するとコードが自動でビルドされます．
build.shのTP_DIRをTextParseをインストールしたパスに変更してください．
CMAKE_INSTALL_PREFIXは実行ファイルが生成させるディレクトリのパスです．
実行する際はすべてのinputファイルが存在するディレクトリ上でTwoDimensionalDiffusionを実行してください．
コマンドライン引数としてtest.tpを与えてください

# output file
- out_Cディレクトリに全ての相の結果を足し合わせた濃度変化ファイル(vtu)が生成されます．
- fluidディレクトリに液体領域の濃度変化ファイル(vtu)が生成されます．
- solidディレクトリに固体領域の濃度変化ファイル(vtu)が生成されます．
- vesselディレクトリに血管領域の濃度変化ファイル(vtu)が生成されます．
- fluid_c.h5 : 液体領域の濃度の経時変化情報が格納されています．
- solid_c.h5 : 固体領域の濃度の経時変化情報が格納されています．
- vessel_c.h5 : 血管領域の濃度の経時変化情報が格納されています．
- sum_O17.h5 : 全ての相の濃度の経時変化情報が格納されています．