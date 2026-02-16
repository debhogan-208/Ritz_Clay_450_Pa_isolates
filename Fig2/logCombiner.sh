# Kenneth B. Hoehn
# 2/10/26
# Re-sample and annotate BEAST results


../beast/bin/logcombiner -renumber -resample 100000 -log sample_full.log backup/sample_full.log -o sample_full_combined_100000.log

../beast/bin/logcombiner -renumber -resample 100000 -log lobe_tree_with_trait.trees backup/lobe_tree_with_trait.trees -o lobe_tree_with_trait_combined_100000.trees

../beast/bin/logcombiner -renumber -resample 100000 -log sample_full.trees backup/sample_full.trees -o sample_full_combined_100000.trees

../beast/bin/logcombiner -renumber -resample 100000 -log tree_asr.trees backup/tree_asr.trees -o tree_asr_combined_100000.trees

../beast/bin/treeannotator -b 0 -lowMem  tree_asr_combined_100000.trees tree_asr_combined_100000.tree

../beast/bin/treeannotator -b 0 -lowMem  lobe_tree_with_trait_combined_100000.trees lobe_tree_with_trait_combined_100000.tree

../beast/bin/treeannotator -b 0 -lowMem  sample_full_combined_100000.trees sample_full_combined_100000.tree

