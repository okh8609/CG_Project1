# CG_Project1
NTUST Computer Graphics Project1 Image Processing
## 操作介紹
![](https://i.imgur.com/c9H2DrS.png)
## 相關技術
我有實作以下的技術：
#### 轉灰階
RGB給不同的權重，產出灰階值。
#### Uniform Quantization
平均去取樣。
#### Populosity Quantization
以最受歡迎的幾種顏色做代表，這個我有用OpenMP去加速。
#### Naive Dithering
二值化。
#### Brightness Dithering
根據平均分布做二值化。
#### Random Dithering
增加雜訊，在做二值化。
#### Cluster Dithering
根據mask做二值化。
#### Floyd Dithering
把自己的誤差推往旁邊的二值化。
#### Color Floyd Dithering
RGB分別做Floyd Dithering。
#### Box Filter
跟旁邊的像素做平均。
#### Barlette Filter
帶權重的跟旁邊像素做平均，越遠權重越小，呈現性關係。
#### Gaussian Filter
帶權重的跟旁邊像素做平均，越遠權重越小，呈高斯函數分布。
#### Edge Detect Filter
去除低頻的資訊，就只留下高頻的資訊。
#### Edge Enhance Filter
不斷地強調高頻的資訊。
#### Half Size
把圖片大小便一半。
#### Double Size
把圖片變兩倍大。

## 畫面呈現

![](https://i.imgur.com/V9EGrDi.png)

![](https://i.imgur.com/ZfjYf5s.png)

![](https://i.imgur.com/gRy8Vih.png)
