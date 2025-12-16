# Orbital Mechanics Library

軌道力学、座標変換、摂動モデルなどに関する包括的なライブラリです。

## 主な機能
- **軌道要素変換**: `oe2rv`, `rv2oe`, `mee2coe` など
- **座標変換**: `itrf2gcrf`, `teme2J2000`, `ecef2LatLonH` など (IAU2006/2000A 歳差・章動モデル等を含む)
- **重力モデル**: `egm2008`, `earthG` (球調和関数展開)
- **その他**: ケプラー方程式を解く `keplerEq`、日食判定 `shadow` など
