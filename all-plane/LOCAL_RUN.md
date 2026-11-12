# ローカルでのビルドと検証

対象はmacOSまたはLinuxです。Python 3.9以上と、128ビット整数に対応するClangまたはGCCが必要です。Pythonの追加パッケージはありません。スクリプトはソフトウェアを自動インストールしません。

本パッケージを展開した `rational_square_search` ディレクトリで実行します。

```sh
python3 local_run.py check
python3 local_run.py build
python3 local_run.py smoke
```

`check` は配布ファイルを `SHA256SUMS` と照合します。`build` は `CXX` があればその指定を使い、なければ `clang++`、`g++`、`c++` の順で探します。元の探索ソースを変更せず、ソースと同じディレクトリに `all_plane_gap_v2_local`、`all_plane_gap_v3_local`、`all_plane_gap_v4_local`、`square_search_local` の4本を作成します。

`smoke` は既定のv4で全平面の分母10万の生成器の内部検証を実行し、2ワーカーの分割実行も検査します。診断用の3距離候補は、v4専用の参照CSVの非空集合と比較します。対応する参照JSONに記録された厳密再構成組数・診断点数・4距離解数とも照合します。v4の追加必要条件は以前の3距離候補を正当に除くことがあるため、旧版の10点を保存することは要求しません。

凍結した参照は `research_v4/final_100000.json` と `research_v4/final_three_100000.csv` です。分母10万の厳密再構成組数は **8,710**、診断用3距離点は **7点**、4距離解は **0件**です。参照ファイルのハッシュもsmokeの記録に保存します。

さらに内部版の分母10万の結果を同梱の記録と比較します。時間・メモリの計測値は機種により変わるので一致条件に含めません。分母200万や内部1億の大きな探索は自動起動しません。

## 小さな全平面探索を実行する

```sh
python3 local_run.py run --bound 100000 --workers 2
```

既定の探索器はv4です。旧版との比較には `--engine v2` または `--engine v3` を指定できます。各版の研究結果と適用範囲は [README.md](README.md) を参照してください。

結果は毎回新しい `local_runs/run_日時_識別子/` に作成します。過去の結果や同梱のログを上書きしません。完走時にはその中の `shards/summary.json` の場所を表示します。`local_build.json` は現在のローカルビルドの記録で、過去のビルド記録は `local_runs/build_日時_識別子/` に残ります。

この起動スクリプトは最大2ワーカー、1回50秒の実時間期限、子プロセスごと30CPU秒、1ファイル16MiBの制限を設定します。Linuxでは子プロセスごと512MiBの仮想メモリ制限も設定します。macOSでは仮想メモリ量と実使用RAMが大きく異なるため、その制限は設定しません。既定の小規模検証は以前のローカル検証予算（120CPU秒、2GiB RAM、1GiB出力）に収まるように選んでいます。

`--bound` の数値上限は200万ですが、計算時間の制限に達した場合は停止し、その範囲を「完走」や「解なし」と報告しません。200万までの完了済み記録は `sweep_v2_2000000/` に入っています。これは旧v2の結果であり、今回のv4の小規模検証を200万までの再検証と扱ってはいけません。既存の大規模結果の再計算は配備確認には不要です。

v4は `all_plane_gap_v4.cpp` から `all_plane_gap_v3.cpp`、その先の `all_plane_sieved.cpp` を参照し、`research_v4_local/lifted_projection.hpp` と `research_v4_math/critical_extension.hpp` も使います。ヘッダーを含めた全依存ファイルと実行ドライバーのハッシュを保持します。単独のcppファイルだけをコピーしてビルドしないでください。

`run_sharded_v4.py` は既定の全追加条件が有効であることも検査します。比較用に条件を無効化した探索器や旧版を渡した場合、v4完走の記録として受け入れません。検証が省略されないよう、Pythonの `-O` と `PYTHONOPTIMIZE` が有効な単独起動も拒否します。

## macOSについて確認した点

- 全平面版はC++17、内部版はC++20でビルドします。`__int128` を使うためClang/GCC系のコンパイラが必要です。
- `run_sharded_v2.py`、`run_sharded_v3.py`、`run_sharded_v4.py` の `resource` はPOSIX向けです。並列処理はOSの子プロセスを使い、OpenMPなどの追加ライブラリには依存しません。
- 平方根の浮動小数計算は初期値だけです。補正と最終的な平方判定は整数演算です。
- 内部版 `search.cpp` の `peak_rss_kib` はLinuxの単位を前提としたフィールド名です。macOSの `getrusage().ru_maxrss` はバイトなので、その値をKiBとして比較してはいけません。この違いは距離の判定や候補数には影響しません。既存のソース・ハッシュを保つため、ここではソースを変更していません。

このパッケージを作った環境でのLinux検証と、ユーザーのMacでの実行は別です。ユーザーのMacで `check`・`build`・`smoke` が成功した受領記録が得られるまでは、Macへの配備・検証が完了したとは扱いません。`local_runs/smoke_日時_識別子/smoke.json` に、実行したOS・CPUアーキテクチャ・ソースハッシュと検査結果を記録します。

研究の到達点・未実装の改善候補は `README.md`、`literature_review_2026-11-04.md`、`prior_arithmetic_overlap.md` を参照してください。
