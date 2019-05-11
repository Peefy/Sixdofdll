# Sixdofdll
 六自由度平台核心算法库

包括六自由度平台正解算法、反解算法、洗出算法，可以灵活设置平台机械尺寸参数和算法参数

## Api函数

函数名称|函数功能
-|-
Control|反解算法：由姿态角获取缸的伸长行程，位移单位为毫米(mm)，角度单位为角度(deg),
FromLengthToPose|正解算法：由缸的伸长行程获取姿态角，位移单位为毫米(mm)，角度单位为角度(deg)
GetTopPosition|获取顶部铰链的空间坐标，坐标原点为平台升起后中立位顶部平台中心，单位为毫米(mm)
GetBottomPosition|获取底部铰链的空间坐标，坐标原点为平台升起后中立位顶部平台中心，单位为毫米(mm)
SetTopPosition|设置顶部铰链的空间坐标，坐标原点为平台升起后中立位顶部平台中心，单位为毫米(mm)
SetBottomPosition|设置底部铰链的空间坐标，坐标原点为平台升起后中立位顶部平台中心，单位为毫米(mm)
SetPlatformPara|设置六自由度平台的结构参数，单位为毫米(mm)
WashOutFiltering|六自由度平台的洗出算法
SetWashOutFilterPara|设置六自由度平台的洗出算法的参数

## Thanks

[cminpack](https://github.com/devernay/cminpack)

