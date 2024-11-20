I have this code.

print("---------- 🥰🥰🥰🥰 RUN MERF 🥰🥰🥰🥰 ----------")
mrf = MERF()
print("X columns: ", X.columns.to_list())
print("Z: ", Z)
print("Y: ", Y)
print("clusters_train: ", clusters_train)

mrf.fit(X, Z, pd.Series(clusters_train), Y)

When I run it I get this.

---------- 🥰🥰🥰🥰 RUN MERF 🥰🥰🥰🥰 ----------
X columns:  ['Unnamed: 0_tax', 'g__Parabacteroides_B_862066', 'g__Coprenecus', 'g__Butyricimonas', 'g__Odoribacter_865974', 'g__Alistipes_A_871404', 'g__Paramuribaculum', 'g__Alistipes_A_871400', 'g__Barnesiella', 'g__Coprobacter', 'g__Phocaeicola_A_858004', 'g__Bacteroides_H', 'g__Prevotella', 'g__Paraprevotella', 'g__Methanobrevibacter_A', 'g__DTU012', 'g__Escherichia_710834', 'g__Parasutterella', 'g__Sutterella', 'g__Haemophilus_D_735815', 'g__Enterobacter_B_713587', 'g__Akkermansia', 'g__Eubacterium_O_258270', 'g__Anaerofustis', 'g__Peptococcus', 'g__QAMH01', 'g__Senegalimassilia', 'g__Adlercreutzia_404257', 'g__Slackia_A', 'g__Eggerthella', 'g__CAG.1427', 'g__Gordonibacter', 'g__Collinsella', 'g__Holdemania', 'g__Longibaculum', 'g__Catenibacterium', 'g__Erysipelatoclostridium', 'g__Faecalibacillus', 'g___2', 'g__Holdemanella', 'g__Merdibacter', 'g__Clostridium_AQ', 'g__Amedibacillus', 'g__Longicatena', 'g__Dielma', 'g__Pauljensenia', 'g__Bifidobacterium_388775', 'g__Acidaminococcus', 'g__Phascolarctobacterium_A', 'g__Veillonella_A', 'g__Dialister', 'g__Streptococcus', 'g__Lactococcus_A_346120', 'g__Lacticaseibacillus', 'g__Turicibacter', 'g__Gemella', 'g__Limiplasma', 'g__SFMI01', 'g___3', 'g__Onthenecus', 'g__Copromorpha', 'g__BX12', 'g___4', 'g__Romboutsia_B', 'g__CCUG.7971', 'g__Intestinibacter', 'g__Terrisporobacter', 'g___5', 'g__PeH17', 'g__Clostridium_T', 'g__CAG.269', 'g___6', 'g__CAG.273', 'g__Merdicola', 'g__CAG.41', 'g__Monoglobus', 'g__ER4', 'g___7', 'g__Dysosmobacter', 'g__Vescimonas', 'g__CAG.83', 'g__Limivicinus', 'g__Faecousia', 'g__Agathobaculum', 'g__Lawsonibacter', 'g__Intestinimonas', 'g__Pseudobutyricicoccus', 'g__Evtepia', 'g__Onthomonas', 'g___8', 'g__RUG13077', 'g__Avispirillum', 'g__SFLA01', 'g__WQUU01', 'g__UMGS1071', 'g__Clostridium_A', 'g__Acutalibacter', 'g__UBA1417', 'g__Anaeromassilibacillus', 'g__Hydrogeniiclostridium', 'g__Ruminococcus_E', 'g__UBA5905', 'g__Bittarella', 'g__Fimenecus', 'g__Faecalibacterium', 'g__Gemmiger_A_73129', 'g__Ruthenibacterium', 'g__Ruminiclostridium_E', 'g___9', 'g__CAG.177', 'g__CAG.217', 'g__Angelakisella', 'g__Massilioclostridium', 'g__Phocea', 'g__Anaerotruncus', 'g__Negativibacillus', 'g__UBA1394', 'g__Ruminococcus_C_58660', 'g__Ruminococcus_D', 'g___10', 'g___11', 'g__Anaerotignum_189125', 'g__CAG.274', 'g__Anaerobutyricum', 'g__Eubacterium_I', 'g__Roseburia', 'g__Agathobacter_164117', 'g__Butyribacter', 'g__RUG115', 'g__CAG.45', 'g__Eubacterium_J', 'g__Eubacterium_G', 'g__Howardella', 'g___1', 'g__Anaerostipes', 'g__Coprococcus_A_121497', 'g__CAG.127', 'g__Lachnospira', 'g__Eubacterium_F', 'g__Coprococcus_A_187866', 'g__Catenibacillus', 'g__Frisingicoccus', 'g__Merdisoma', 'g__Marvinbryantia', 'g__Blautia_A_141780', 'g__Blautia_A_141781', 'g__UMGS1375', 'g__Limivivens', 'g__UBA9414', 'g__Oliverpabstia', 'g__Mediterraneibacter_A_155507', 'g__COE1', 'g__Acetatifactor', 'g__CAG.95', 'g__Hungatella_A_128155', 'g__Hungatella_A_127239', 'g__Sellimonas', 'g__Eisenbergiella', 'g__UBA3402', 'g__Scatomonas', 'g__Clostridium_Q_134516', 'g__Faecalimonas', 'g__Ruminococcus_B', 'g__Clostridium_AP', 'g__CAG.317_146760', 'g__Bariatricus', 'g__Muricomes_149725', 'g__Lachnoclostridium_B', 'g__Schaedlerella', 'g__Mediterraneibacter_A_155590', 'g__Lactonifactor', 'g__Enterocloster', 'g__Ventrisoma', 'g__Clostridium_Q_135853', 'g__Clostridium_Q_135822', 'g__Porcincola', 'g__Copromonas', 'g__Ventrimonas', 'g__Dorea_A', 'g__Massilistercora', 'x_t', 'Unnamed: 0_long']
Z:  [[1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]
 [1.]]
Y:  [32.244781   34.11720866 35.42716527 30.65216963 28.59405902 28.2562908
 33.53394204 39.44467518 43.352564   34.16130722 29.40162853 29.51156263
 37.08516893 27.18705049 29.98377268 31.91106484 31.88517378 27.85409926
 33.24970996 39.8711175  38.24032888 30.26298124 34.06215479 28.08931035
 29.71224227 33.67518842 28.37786514 29.94801741 38.62990675 35.37812846
 44.16388985 36.8178737  33.97265981 32.33288034 32.87667481 38.61588109
 32.12637387 29.18778802 37.61665531 31.0498484  28.94347481 29.88557903
 31.46756393 33.12418673 37.15487084 34.98886855 36.42649257 37.24904506
 37.57444482 38.24905053 38.29771813 26.80131708 33.16394465 37.85531782
 33.98791001 37.09584134 28.02243158 38.5272789  39.04013252 35.01318998
 41.80984656 31.78832715 29.71074294 34.05607592 37.77969107 33.79769745
 35.15347901 46.43449077 32.31623964 35.06202875 34.45716044 35.2005934
 41.80143871 33.98247831 39.67941194 41.78635866 27.29200763 29.27480756
 37.59393845 28.63932413 38.0829356  33.64209085 37.83373797 36.05274929
 30.21951142 40.3703385  33.23102313 37.1765976  29.36985044 38.47984471
 32.11531212 31.84461234 30.69113943 28.64582805 34.33587848 34.69726497
 36.90907287 43.88602846 26.77663542 27.1409109  39.8357322  29.76295347
 37.10876128 32.56626272 29.90705531 31.05556243 27.357311   34.24973973
 31.5975537  36.56626658 34.60207514 34.89554816 31.23990811 27.28386116
 30.14131423 30.14579725]
clusters_train:  [144  53 105  71  55   9  11 139 163 128 157  13 129  52  83 167 193  14
 201 104 156 164 174  57 189  26  73  77  49  66 161 130 177  92 135  45
 182 184  98 188  10 109 140  27  93 127  22 116  32 126 162  85  34 199
  29 170 196  42  72 122 141  47  15 148 147   3 101 111 198  17 142 119
 187  88  65 107  68 118  76 117  91 112 143 133 152  48  69 171  67 149
  21 166 114  23 132   6  20  43 150  80  62 195 178  28 115  86  16 172
 113 146  56  41   5 103  97  96]
---------------------------------------------------------------------------
ValueError                                Traceback (most recent call last)
Cell In[29], line 8
      5 print("Y: ", Y)
      6 print("clusters_train: ", clusters_train)
----> 8 mrf.fit(X, Z, pd.Series(clusters_train), Y)

File ~/anaconda3/lib/python3.10/site-packages/merf/merf.py:219, in MERF.fit(self, X, Z, clusters, y, X_val, Z_val, clusters_val, y_val)
    216 assert len(y_star.shape) == 1
    218 # Do the fixed effects regression with all the fixed effects features
--> 219 self.fe_model.fit(X, y_star)
    220 f_hat = self.fe_model.predict(X)
    222 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ M-step ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

File ~/anaconda3/lib/python3.10/site-packages/sklearn/ensemble/_forest.py:345, in BaseForest.fit(self, X, y, sample_weight)
    343 if issparse(y):
    344     raise ValueError("sparse multilabel-indicator for y is not supported.")
--> 345 X, y = self._validate_data(
    346     X, y, multi_output=True, accept_sparse="csc", dtype=DTYPE
    347 )
    348 if sample_weight is not None:
    349     sample_weight = _check_sample_weight(sample_weight, X)

File ~/anaconda3/lib/python3.10/site-packages/sklearn/base.py:565, in BaseEstimator._validate_data(self, X, y, reset, validate_separately, **check_params)
    563         y = check_array(y, input_name="y", **check_y_params)
    564     else:
--> 565         X, y = check_X_y(X, y, **check_params)
    566     out = X, y
    568 if not no_val_X and check_params.get("ensure_2d", True):

File ~/anaconda3/lib/python3.10/site-packages/sklearn/utils/validation.py:1106, in check_X_y(X, y, accept_sparse, accept_large_sparse, dtype, order, copy, force_all_finite, ensure_2d, allow_nd, multi_output, ensure_min_samples, ensure_min_features, y_numeric, estimator)
   1101         estimator_name = _check_estimator_name(estimator)
   1102     raise ValueError(
   1103         f"{estimator_name} requires y to be passed, but the target y is None"
   1104     )
-> 1106 X = check_array(
   1107     X,
   1108     accept_sparse=accept_sparse,
   1109     accept_large_sparse=accept_large_sparse,
   1110     dtype=dtype,
   1111     order=order,
   1112     copy=copy,
   1113     force_all_finite=force_all_finite,
   1114     ensure_2d=ensure_2d,
   1115     allow_nd=allow_nd,
   1116     ensure_min_samples=ensure_min_samples,
   1117     ensure_min_features=ensure_min_features,
   1118     estimator=estimator,
   1119     input_name="X",
   1120 )
   1122 y = _check_y(y, multi_output=multi_output, y_numeric=y_numeric, estimator=estimator)
   1124 check_consistent_length(X, y)

File ~/anaconda3/lib/python3.10/site-packages/sklearn/utils/validation.py:879, in check_array(array, accept_sparse, accept_large_sparse, dtype, order, copy, force_all_finite, ensure_2d, allow_nd, ensure_min_samples, ensure_min_features, estimator, input_name)
    877         array = xp.astype(array, dtype, copy=False)
    878     else:
--> 879         array = _asarray_with_order(array, order=order, dtype=dtype, xp=xp)
    880 except ComplexWarning as complex_warning:
    881     raise ValueError(
    882         "Complex data not supported\n{}\n".format(array)
    883     ) from complex_warning

File ~/anaconda3/lib/python3.10/site-packages/sklearn/utils/_array_api.py:185, in _asarray_with_order(array, dtype, order, copy, xp)
    182     xp, _ = get_namespace(array)
    183 if xp.__name__ in {"numpy", "numpy.array_api"}:
    184     # Use NumPy API to support order
--> 185     array = numpy.asarray(array, order=order, dtype=dtype)
    186     return xp.asarray(array, copy=copy)
    187 else:

File ~/anaconda3/lib/python3.10/site-packages/pandas/core/generic.py:2070, in NDFrame.__array__(self, dtype)
   2069 def __array__(self, dtype: npt.DTypeLike | None = None) -> np.ndarray:
-> 2070     return np.asarray(self._values, dtype=dtype)

ValueError: could not convert string to float: 'AAL-144.BL'


​