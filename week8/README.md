For Step 1.1: Rscript runChicago.R ./raw/PCHIC_Data/GM_rep1.chinput,./raw/PCHIC_Data/GM_rep2.chinput,./raw/PCHIC_Data/GM_rep3.chinput ./output --design-dir ./raw/Design/ --en-feat-list ./raw/Features/featuresGM.txt --export-format washU_text

Step 1.2: CTCF is enriched, which is expected for active promoters and enhancers. This makes sense if the bait sequence is an active promoter. H3K4me1 is enriched, which is expected for active or primed enhancers. This interaction makes sense if the bait sequence is an active promoter. H3K4me3 is enriched, this is expected at the transcription start site of active genes. This makes sense if the bait sequence is an active promoter. H3k27ac is enriched, which is expected at the transcription start site of active genes. This makes sense if the bait sequence is an active promoter. H3K27me3 is not enriched, this normally indicates gene repression. This makes sense if the bait sequence is an active promoter. H3K9me3 is enriched, which is associated with heterochromatin. This is unexpected as the other enriched markers indicate that the bait sequence is an active promoter.

For Step 2.2

Top 6 promoter-promoter interactions
1. chr20 44438565 44565593 . 1000 34.77 . 0 chr20 44438565 44442365                                           UBE2C + chr20 44562442 44565593 PCIF1 +
2. chr20 44438565 44607204 .  986 34.29 . 0 chr20 44438565 44442365                                           UBE2C + chr20 44596299 44607204 FTLP1;ZNF335 +
3. chr21 26837918 26939577 .  978 34.02 . 0 chr21 26837918 26842640                                          snoU13 + chr21 26926437 26939577 MIR155HG +
4. chr20 44452862 44565593 .  974 33.89 . 0 chr20 44452862 44471524                                     SNX21;TNNC2 + chr20 44562442 44565593 PCIF1 +
5. chr20 17660712 17951709 .  973 33.85 . 0 chr20 17660712 17672229                                           RRBP1 + chr20 17946510 17951709 MGME1;SNX5 +
6. chr20 24972345 25043735 .  973 33.84 . 0 chr20 24972345 24985047                                           APMAP + chr20 25036380 25043735 ACSS1 +

Top 6 promoter-enhancer interactions
1. chr21 26797667 26939577 .  952 33.13 . 0 chr21 26926437 26939577                                        MIR155HG + chr21 26797667 26799364 . -
2. chr20 55957140 56074932 .  928 32.29 . 0 chr20 55957140 55973022                              RBM38;RP4-800J21.3 + chr20 56067414 56074932 . -
3. chr21 26790966 26939577 .  838 29.17 . 0 chr21 26926437 26939577                                        MIR155HG + chr21 26790966 26793953 . -
4. chr20  5585992  5628028 .  830 28.88 . 0 chr20  5585992  5601172                                          GPCPD1 + chr20  5625693  5628028 . -
5. chr21 26793954 26939577 .  754 26.23 . 0 chr21 26926437 26939577                                        MIR155HG + chr21 26793954 26795680 . -
6. chr20  5515866  5933156 .  750 26.08 . 0 chr20  5929472  5933156                                      MCM8;TRMT6 + chr20  5515866  5523933 . -

For Step 2.3
	It makes sense for RBM38 to be interacting with enhancers as it is normally highly expressed in bone marrow and lymph nodes, indicating that it is expressed in immune cells which GM12878 is derived from.
	It makes sense for MIR155HG to be interacting with enhancers as it is often expressed at high levels in lymphoma and GM12878 was immortalized using EBV which is a common cause of lymphoma.