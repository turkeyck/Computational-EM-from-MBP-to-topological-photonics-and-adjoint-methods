# 計算電磁學與拓撲光子學：從平面波展開到可微分 RCWA

Oct 29, 2026 · @jr lu

## 課程概述

本講義用 11 講把兩年計畫的技術內容串成一門研究所課程：從光子能帶的本徵問題出發，算出拓撲不變量，再走到 RCWA、FDTD／GPU 與可微分逆向設計。水準對標史丹佛 EE 研究所的計算電磁學與光子學課程，每講都以「推導 → 演算法 → 可驗證的程式作業」收尾。

**先修**：電磁學（研究所程度）、線性代數（Hermitian 算子、廣義本徵問題）、傅立葉分析、Python／NumPy。量子力學只需大學程度；第 5 講的 QFT 為伴讀，不影響主線。

**設計原則**

- 每一個數值結果都要有獨立的對照答案：解析解、MPB、S4，或能量守恆。
- 先求正確，再求快：Python 原型 → 測試 → C++／GPU 改寫。
- 每講的作業都是公開 repo 的一個可執行範例，同時也是作品集。

| 講次 | 主題 | 對應階段 | 建議時程 | 交付物 |
| --- | --- | --- | --- | --- |
| 1–2 | 本徵問題與 PWE | 第一階段 | 第 1–5 週 | 方晶格能帶，對照 MPB |
| 3–4 | Berry 曲率、Chern 數、邊緣態 | 第一階段 | 第 6–12 週 | 非零 Chern 數 + 超晶胞邊緣態 |
| 5 | 量子異常（伴讀） | 第一階段起並行 | 全程 | 讀書筆記 |
| 6–8 | RCWA、Li 法則、S 矩陣 | 第二階段 | 第 4–7 個月 | 1D/2D RCWA，對照 S4 |
| 9 | C++ 核心與 pybind11 | 第二階段 | 第 7–9 個月 | pip 可安裝套件 |
| 10 | FDTD 與 GPU | 第二階段（公司工作） | 並行 | 加速倍率報告（不含機密） |
| 11 | 可微分求解器與 adjoint | 第三階段 | 第 9–18 個月 | 逆向設計示範 |

**每週節奏（8–10 小時）**：平日 3 個晚上各 1.5 小時讀推導與寫程式，週末一個 4 小時區塊做整合、驗證與寫筆記。

**評量**：每講的「驗收標準」就是期末考。達不到就不往下走，回頭縮小範圍。

## 第 1 講：馬克士威方程作為 Hermitian 本徵問題

核心結論：在線性、無損、無色散的介質中，光子晶體的模態就是一個 Hermitian 算子的本徵態，所以頻率是實數、模態正交，能帶理論的整套工具都能直接搬過來。

**學習目標**：寫出主方程並證明其 Hermitian 性；理解 Bloch 定理如何把無限問題化成單一晶胞；能用尺度不變性換算設計。

### 1.1 主方程

假設時諧場 e^(−iωt)、無自由電荷與電流、μ = μ₀、ε(r) 為實數。消去 E 得到磁場的主方程：

```latex
\Theta \mathbf{H} \equiv \nabla \times \left( \frac{1}{\varepsilon(\mathbf{r})} \nabla \times \mathbf{H} \right) = \left( \frac{\omega}{c} \right)^2 \mathbf{H}, \qquad \nabla \cdot \mathbf{H} = 0
```

在內積 ⟨F, G⟩ = ∫ F\*·G d³r 下，Θ 是 Hermitian 且半正定。兩個直接推論：ω² 為非負實數；不同頻率的模態彼此正交。橫向條件 ∇·H = 0 不能省：它排除了 ω = 0 的縱向假解，數值上也會以偽模態的形式回來找你。

### 1.2 Bloch 定理與布里淵區

若 ε(r + R) = ε(r) 對所有晶格向量 R 成立，Θ 與平移算子對易，模態可寫成

```latex
\mathbf{H}_{n\mathbf{k}}(\mathbf{r}) = e^{i\mathbf{k}\cdot\mathbf{r}}\, \mathbf{u}_{n\mathbf{k}}(\mathbf{r}), \qquad \Theta_{\mathbf{k}} = (\nabla + i\mathbf{k}) \times \frac{1}{\varepsilon} (\nabla + i\mathbf{k}) \times
```

u 只需在一個晶胞內求解，k 只需取遍第一布里淵區；能帶 ω\_n(k) 就是這族本徵值隨 k 的變化。第 3 講的 Berry 曲率，正是定義在 u\_nk 隨 k 變化的幾何上。

### 1.3 二維結構：TE 與 TM 分離

對 z 方向均勻的 2D 晶體且 k\_z = 0，場分成兩個獨立的純量問題。本講義採 Joannopoulos 的慣例：TM 指電場沿 z。

```latex
\text{TM } (E_z):\; -\nabla^2 E_z = \left(\tfrac{\omega}{c}\right)^2 \varepsilon(\mathbf{r})\, E_z \qquad\qquad \text{TE } (H_z):\; -\nabla \cdot \left( \tfrac{1}{\varepsilon} \nabla H_z \right) = \left(\tfrac{\omega}{c}\right)^2 H_z
```

TM 是廣義 Hermitian 本徵問題（右邊有 ε 權重），TE 是標準 Hermitian 問題但算子裡有 1/ε。這個差異在第 2 講決定了傅立葉展開要怎麼做。

### 1.4 尺度不變性與變分原理

若 ε′(r) = ε(r/s)，則 ω′ = ω/s：沒有基本長度尺度，所以能帶一律用無因次頻率 ωa/2πc 表示。變分原理說明最低能帶會把電場能量集中在高 ε 區域，這是判讀能隙出現在哪一個偏振、哪一對能帶之間的直覺來源。

### 習題

1. 證明 Θ 在上述內積下為 Hermitian；指出哪一步用到 ε 為實數，並說明有損耗時會壞在哪裡。
2. 從 Θ 推導 Θ\_k，並證明 ω\_n(k) = ω\_n(k + G)。
3. 一維 Bragg 堆疊（折射率 n₁、n₂，厚度 d₁、d₂，週期 Λ）正入射時，用轉移矩陣推出 Rytov 色散關係：

```latex
\cos(k\Lambda) = \cos(k_1 d_1)\cos(k_2 d_2) - \frac{1}{2}\left( \frac{n_1}{n_2} + \frac{n_2}{n_1} \right) \sin(k_1 d_1)\sin(k_2 d_2)
```

4. 對四分之一波堆疊，求第一能隙的中心頻率與相對寬度，並解釋為何 n₁/n₂ → 1 時能隙趨近於零。

**驗收標準**：習題 3 的解析能帶會成為第 2 講 PWE 程式的第一個單元測試。

**閱讀**：Joannopoulos 等人《Photonic Crystals: Molding the Flow of Light》第 2 版（2008）第 2–3 章；Sakoda《Optical Properties of Photonic Crystals》（Springer, 2001）。

## 第 2 講：平面波展開法（PWE）

核心結論：把場與介電函數都展開成倒晶格平面波，主方程就變成一個有限維的（廣義）Hermitian 矩陣本徵問題；精度的關鍵不在平面波數量，而在介電函數的傅立葉係數怎麼取。

**學習目標**：推導 TM／TE 的 PWE 矩陣；理解 Laurent 規則與「逆規則」的收斂差異；寫出可畫 Γ–X–M–Γ 能帶圖的程式。

### 2.1 展開與矩陣形式

取 N 個倒晶格向量 G（截斷在 |G| < G\_max），令

```latex
E_z(\mathbf{r}) = \sum_{\mathbf{G}} e_{\mathbf{G}}\, e^{i(\mathbf{k}+\mathbf{G})\cdot\mathbf{r}}, \qquad \varepsilon(\mathbf{r}) = \sum_{\mathbf{G}} \varepsilon_{\mathbf{G}}\, e^{i\mathbf{G}\cdot\mathbf{r}}
```

代入 TM 方程並對 e^(i(k+G)·r) 投影：

```latex
|\mathbf{k}+\mathbf{G}|^2\, e_{\mathbf{G}} = \left(\frac{\omega}{c}\right)^2 \sum_{\mathbf{G}'} \varepsilon_{\mathbf{G}-\mathbf{G}'}\, e_{\mathbf{G}'} \quad\Longleftrightarrow\quad \mathbf{K}\,\mathbf{e} = \lambda\, [\![\varepsilon]\!]\,\mathbf{e}
```

K 是對角矩陣，⟦ε⟧ 是由 ε\_(G−G′) 組成的 Toeplitz 矩陣（Hermitian 正定）。直接用 scipy.linalg.eigh(K, Eps) 解廣義問題，保留 Hermitian 結構，不要先求反矩陣再丟給一般的 eig。

TE 方程的投影為

```latex
\sum_{\mathbf{G}'} (\mathbf{k}+\mathbf{G})\cdot(\mathbf{k}+\mathbf{G}')\, \eta_{\mathbf{G}\mathbf{G}'}\, h_{\mathbf{G}'} = \left(\frac{\omega}{c}\right)^2 h_{\mathbf{G}}
```

這裡 η 代表 1/ε 的矩陣，而它有兩種取法，正是本講的重點。

### 2.2 Laurent 規則與逆規則

- **Laurent 規則**：η = ⟦1/ε⟧，直接取 1/ε 的傅立葉係數。
- **逆規則**：η = ⟦ε⟧⁻¹，先組 ε 的 Toeplitz 矩陣再求反。

兩者在 N → ∞ 時一致，但截斷後不同。早期的 PWE 計算就因介電函數傅立葉級數收斂差而吃過苦頭（Sözüer、Haus、Inguva 1992 有系統性討論）；Li（1996）從「兩個不連續函數相乘」的傅立葉分解規則給出解釋：當乘積在介面處連續、因子各自不連續時，必須用逆規則。這和第 7 講 RCWA 的 Li 法則是同一個定理，PWE 是你第一次親手看到它。

### 2.3 介電函數的傅立葉係數

對半徑 R 的圓柱（介電常數 ε\_a，背景 ε\_b，填充率 f = πR²/A）有解析式：

```latex
\varepsilon_{\mathbf{G}} = \varepsilon_b\,\delta_{\mathbf{G},0} + (\varepsilon_a - \varepsilon_b)\, f \cdot \frac{2 J_1(|\mathbf{G}|R)}{|\mathbf{G}|R} \quad (\mathbf{G} \neq 0 \text{ 時 } \delta \text{ 項為零，} \mathbf{G}=0 \text{ 時為平均值})
```

任意形狀則在細網格上取樣後做 FFT，但要小心 Gibbs 振盪與網格對齊造成的鋸齒。進階做法是次像素平滑（Farjadpour 等人，2006），MPB 就是這樣做的。

### 2.4 演算法與成本

1. 選倒晶格基底與 G 集合（以 |G| 排序截斷，而非方形截斷，以保持旋轉對稱）。
2. 預先算好 ε\_(G−G′) 表，對每個 k 只重組對角的 K。
3. 沿高對稱路徑 Γ–X–M–Γ 取 k，求最低 n\_bands 個本徵值，轉成 ωa/2πc。

稠密解法是 O(N³)，N 約數百時筆電可承受。MPB 改用迭代本徵解法，並以 FFT 在實空間與倒空間之間切換來施加算子，單次迭代 O(N log N)，這是你之後改寫 C++ 時的方向。

### 習題

1. 實作 1D 版本，與第 1 講的 Rytov 解析能帶比對，畫出相對誤差隨 N 的變化。
2. 實作 2D TM：方晶格介電柱 ε = 12、R = 0.2a（MPB 官方教學的範例）。MPB 教學輸出的 TM 能隙位於第 1、2 能帶之間，約從 0.2826 到 0.4193（ωa/2πc），相對寬度約 39%；你的結果要與此對照。
3. 實作 2D TE，分別用 Laurent 規則與逆規則，畫出某一能帶頻率對 N 的收斂曲線（log–log），並估計兩者的收斂階數。
4. 證明：當 ε 只在 x 方向不連續、而 D\_x 連續時，⟦ε⟧⁻¹ 才是 1/ε 乘法的正確截斷表示。

**驗收標準**：2D TM 能帶與 MPB 在前 6 條能帶的差異小於 0.5%。

**閱讀**：Ho、Chan、Soukoulis, PRL 65, 3152 (1990)（早期全向量 PWE 經典）；Sözüer、Haus、Inguva, PRB 45, 13962 (1992)（PWE 收斂問題）；Farjadpour 等人, Opt. Lett. 31, 2972 (2006)（次像素平滑）；Li, JOSA A 13, 1870 (1996)；Johnson & Joannopoulos, Optics Express 8, 173 (2001)（MPB 的演算法論文）。

## 第 3 講：Berry 相位、Berry 曲率與 Chern 數

核心結論：Chern 數是 Berry 曲率在整個布里淵區上的積分，必為整數；數值上用 Fukui–Hatsugai–Suzuki（FHS）的格點方法，可以在相當粗的 k 網格上得到「精確的」整數，而且完全不受本徵解法給出的任意相位影響。

**學習目標**：定義 Berry 聯絡與曲率並理解其規範自由度；推導 FHS 方法；先在已知答案的兩能帶模型上驗證，再接到 PWE。

### 3.1 定義

對第 n 條能帶的週期部分 u\_nk，定義 Berry 聯絡、曲率與 Chern 數：

```latex
\mathbf{A}_n(\mathbf{k}) = i\,\langle u_{n\mathbf{k}} | \nabla_{\mathbf{k}} u_{n\mathbf{k}} \rangle, \qquad \Omega_n(\mathbf{k}) = \partial_{k_x} A_{n,y} - \partial_{k_y} A_{n,x}, \qquad C_n = \frac{1}{2\pi} \int_{\mathrm{BZ}} \Omega_n \, d^2k \in \mathbb{Z}
```

光子問題的內積必須和本徵問題一致：TM（E\_z）形式是廣義問題，內積要帶 ε 權重，⟨u|v⟩ = ∫ ε u\* v；用 H 場形式則是一般內積。內積取錯，曲率的數值會錯，但積分出來通常仍是整數，所以這個錯很容易漏掉。

規範變換 u → e^(iφ(k)) u 讓 A 改變 ∇φ，但 Ω 不變；布里淵區是環面（沒有邊界），所以 C 是規範不變的拓撲量。整數性的直觀理由：C 衡量的是無法在整個環面上選出連續相位的程度。

### 3.2 對稱性的限制

- 時間反演對稱（實數、互易的 ε 與 μ）：Ω(k) = −Ω(−k)，因此 C = 0。
- 空間反演對稱：Ω(k) = Ω(−k)。
- 兩者同時存在：在非簡併處 Ω ≡ 0。

所以要得到非零 Chern 數，必須打破時間反演，這就是第 4 講要用旋磁材料的原因。只打破空間反演（例如 valley 光子晶體）時，總 Chern 數仍為零，但 Ω 會集中在 K 與 K′ 附近，可以定義 valley Chern 數。

### 3.3 Fukui–Hatsugai–Suzuki 方法

在 L₁ × L₂ 的 k 網格上，對每個格點與方向 μ 定義連結變數（U(1) link）：

```latex
U_\mu(\mathbf{k}) = \frac{\langle u(\mathbf{k}) | u(\mathbf{k}+\hat{\mu}) \rangle}{|\langle u(\mathbf{k}) | u(\mathbf{k}+\hat{\mu}) \rangle|}, \qquad F_{12}(\mathbf{k}) = \operatorname{Im} \ln \left[ U_1(\mathbf{k})\, U_2(\mathbf{k}+\hat{1})\, U_1(\mathbf{k}+\hat{2})^{-1}\, U_2(\mathbf{k})^{-1} \right]
```

其中 ln 取主值 (−π, π\]，C = (1/2π) Σ\_k F₁₂(k)。每個 F₁₂ 都是規範不變的（任意相位在乘積中互相抵消），而總和恆為整數；只要網格夠細、使每個小格子的通量遠小於 π，結果就等於真正的 Chern 數。

**簡併能帶**：若某些能帶在 BZ 內相交，就要把它們當成一組，把內積換成重疊矩陣的行列式：U\_μ = det M / |det M|，其中 M\_mn = ⟨u\_m(k)|u\_n(k+μ̂)⟩。算出來的是整組能帶的總 Chern 數。

### 3.4 PWE 的實作細節：週期規範

環面要閉合，必須讓 u\_(k+G) 與 u\_k 對應同一個態。在平面波基底下，u\_(k+G) = e^(−iG·r) u\_k 等價於係數的指標平移：c\_(G′)(k + G) = c\_(G′+G)(k)。實作時只在網格內部各點解本徵問題；網格邊界上的連結，則用內部點的係數做指標平移得到，而不是在 k + G 處重新求解（重新求解會引入另一個任意相位與不同的截斷集合）。

### 習題

1. 證明 F₁₂ 在 u → e^(iφ) u 下不變，且 Σ F₁₂ / 2π 為整數。
2. 先在 Qi–Wu–Zhang 兩能帶模型 H(k) = sin k\_x σ\_x + sin k\_y σ\_y + (m + cos k\_x + cos k\_y) σ\_z 上實作 FHS。畫出下能帶 Chern 數對 m 的相圖：|m| > 2 時為 0，0 < |m| < 2 時為 ±1，且兩側符號相反。
3. 把 FHS 接到第 2 講的 PWE（TM，方晶格介電柱）。驗證每條孤立能帶的 C = 0，並畫出 Ω(k) 確認 Ω(k) = −Ω(−k)。
4. 把網格從 6 × 6 加到 48 × 48，記錄 C 何時穩定，並解釋為何 FHS 不需要 Ω 本身已收斂。

**驗收標準**：QWZ 相圖正確；時間反演對稱的光子晶體所有孤立能帶都得到 C = 0。

**閱讀**：Fukui、Hatsugai、Suzuki, J. Phys. Soc. Jpn. 74, 1674 (2005)；Xiao、Chang、Niu, Rev. Mod. Phys. 82, 1959 (2010)（Berry 相位綜論）；Asbóth、Oroszlány、Pályi《A Short Course on Topological Insulators》第 2 章（Berry 相位與 Chern 數）與第 6 章（QWZ 模型），arXiv:1509.02295 免費；Qi、Wu、Zhang, PRB 74, 085308 (2006)（QWZ 模型原始論文）。

## 第 4 講：旋磁光子晶體、超晶胞與 bulk–edge 對應

核心結論：加磁場的旋磁材料讓 μ 張量出現虛數的非對角項，打破時間反演，能帶就能帶非零 Chern 數；能隙內單向邊緣模的淨數目，等於能隙兩側塊材 Chern 數總和的差。

**學習目標**：推導含 μ 張量的 TM 方程並放進 PWE；重現 Wang 等人（2008）的 Chern 數；用超晶胞看到邊緣模並驗證 bulk–edge 對應。

### 4.1 旋磁介質的 TM 方程

磁偏置沿 z 的鐵氧體（如 YIG），在微波頻段的磁導率為

```latex
\hat{\mu} = \mu_0 \begin{pmatrix} \mu & i\kappa & 0 \\ -i\kappa & \mu & 0 \\ 0 & 0 & 1 \end{pmatrix}
```

它是 Hermitian 但非實數，所以 μ̂\* ≠ μ̂，時間反演被打破。對 TM 模（E = E\_z ẑ），令 ν = μ̂\_⊥⁻¹（只取面內 2 × 2 部分，以 μ₀ 歸一化），J = \[\[0, 1\], \[−1, 0\]\]，可推得

```latex
-\nabla \cdot \left( \tilde{\nu}\, \nabla E_z \right) = \left( \frac{\omega}{c} \right)^2 \varepsilon(\mathbf{r})\, E_z, \qquad \tilde{\nu} = -J\, \nu\, J
```

各向同性時 ν̃ 退化成 1/μ，回到第 1 講的形式。PWE 矩陣變成

```latex
\sum_{\mathbf{G}'} (\mathbf{k}+\mathbf{G})^{\mathsf{T}}\, \tilde{\nu}_{\mathbf{G}-\mathbf{G}'}\, (\mathbf{k}+\mathbf{G}')\, e_{\mathbf{G}'} = \left(\frac{\omega}{c}\right)^2 \sum_{\mathbf{G}'} \varepsilon_{\mathbf{G}-\mathbf{G}'}\, e_{\mathbf{G}'}
```

左邊是複數 Hermitian 矩陣。ν̃ 的四個分量各自在柱體邊界不連續，第 2 講的 Laurent／逆規則問題在這裡會再出現一次，收斂測試不可省。

### 4.2 基準案例：Wang–Chong–Joannopoulos–Soljačić（2008）

方晶格 YIG 柱置於空氣中，ε = 15、R = 0.11a；在 1600 高斯偏置、4.28 GHz 下，μ = 14、κ = 12.4（皆以真空值歸一化）。不加磁場時，第 2、3 條 TM 能帶在 M 點有二次型簡併；加磁場後簡併打開，形成相對寬度約 10% 的能隙，計入材料色散後約 6%，能帶獲得非零 Chern 數。論文的波導介面是 YIG 晶體與一個轉 45° 的氧化鋁柱方晶格（R = 0.106a、ε = 10），後者所有能帶 Chern 數為零，因此介面上只有一條向前傳播的邊緣模，即能隙以下 Chern 數總和的絕對值為 1。2009 年同組在 Nature 發表微波實驗，觀察到這種抗背向散射的單向傳播。

### 4.3 超晶胞與邊緣模

建一條沿 y 方向 N\_s 個晶胞（例如 20）、沿 x 方向週期的超晶胞，k\_x 掃過 \[−π/a, π/a\]，其餘不變。論文用的是 YIG 與氧化鋁晶體的介面；在 PWE 超晶胞裡，更容易實作的替代方案（本講義的建議，非論文做法）是**磁化反向的疇壁**：上半 κ，下半 −κ。兩側 Chern 數差 2，所以介面上會有兩條同向的手徵模；超晶胞的週期邊界還會在另一端形成第二個疇壁，那裡的兩條模方向相反，判讀時要用場分佈把它們分開。

超晶胞的倒晶格在 y 方向變密（間距 2π/(N\_s a)），平面波數量會隨 N\_s 線性增加，矩陣維度很快上千。這是你第一次感受到稠密 O(N³) 的極限，也是改寫成迭代解法的動機。

### 4.4 Bulk–edge 對應

```latex
N_{\text{right}} - N_{\text{left}} = \sum_{n \le n_{\text{gap}}} C_n^{(\text{A})} - \sum_{n \le n_{\text{gap}}} C_n^{(\text{B})}
```

左邊是介面上向右與向左的邊緣模數目，右邊是兩種塊材在該能隙以下的 Chern 數總和。這條式子就是第 5 講 anomaly inflow 在光子系統的版本：塊材的拓撲決定邊界必須存在無法單獨存在的手徵模。

### 習題

1. 從馬克士威方程推出 4.1 的 TM 方程，並證明 ν̃ 為 Hermitian。
2. 重現 4.2 的能帶與各能帶 Chern 數；與論文的能帶圖對照。
3. 建立疇壁超晶胞，畫出投影能帶，數出能隙內每個疇壁上的邊緣模數與傳播方向，驗證 bulk–edge 對應。
4. （選做）用 Meep（支援旋磁介質的 FDTD）模擬單向邊緣模繞過一個方形缺陷，量測反射率。

**驗收標準**：能隙以下 Chern 數總和的絕對值為 1；疇壁上恰有兩條同向邊緣模。

**閱讀**：Wang 等人, PRL 100, 013905 (2008)；Wang 等人, Nature 461, 772 (2009)；Haldane & Raghu, PRL 100, 013904 (2008)；Lu、Joannopoulos、Soljačić, Nature Photonics 8, 821 (2014)（綜論）。

## 第 5 講（伴讀）：手徵異常、譜流與 anomaly inflow

核心結論：古典守恆的手徵流在量子化後不再守恆，因為沒有任何正規化能同時保住它與規範不變性；單獨一條手徵模的電荷不守恆，必須由高一維塊材的 Chern–Simons 響應「灌入」電流來補償。這正是第 4 講邊緣模存在的場論版本。

本講不影響主線進度，建議與第 1–4 講並行，每週讀 2–3 小時。

### 5.1 1+1 維：譜流圖像

無質量 Dirac 費米子分成右行（E = p）與左行（E = −p）兩支，古典上向量流 j^μ 與軸流 j₅^μ 都守恆。加上均勻電場 E，每個態的動量以 dp/dt = eE 漂移。在長度 L 的系統中，單位動量的態密度為 L/2π，所以 Dirac 海底部的右行態被「抽」上來、左行態被壓下去：

```latex
\frac{d}{dt}\left(N_R - N_L\right) = \frac{eEL}{\pi} \quad\Longrightarrow\quad \partial_\mu j_5^{\mu} = \frac{e}{\pi}\, E = \frac{e}{2\pi}\, \epsilon^{\mu\nu} F_{\mu\nu}
```

總電荷 N\_R + N\_L 守恆，但手徵荷不守恆。關鍵在於無窮深的 Dirac 海：任何截斷（晶格、動量截止）都會把右行與左行連在一起，這就是 Nielsen–Ninomiya 定理的內容。

### 5.2 只有一條手徵模時：anomaly inflow

若邊界上只有右行模，它的電荷本身就不守恆：

```latex
\partial_\mu j^{\mu} = \frac{e}{4\pi}\, \epsilon^{\mu\nu} F_{\mu\nu}
```

這在 1+1 維單獨存在是不一致的。解法是把它放在 2+1 維量子霍爾塊材的邊上：塊材的有效作用量含 Chern–Simons 項，霍爾電流 σxy = C e²/h 垂直流向邊界，恰好補上邊緣的電荷變化（Callan–Harvey, 1985）。

對照第 4 講：光子晶體的 Chern 數 C 對應這裡的 C，疇壁上的單向邊緣模對應這裡的手徵費米子。光子不帶電、沒有真正的「電荷不守恆」，但決定邊緣模數目的拓撲結構是同一個。

### 5.3 3+1 維：ABJ 異常與 Fujikawa 方法

```latex
\partial_\mu j_5^{\mu} = \frac{e^2}{16\pi^2}\, \epsilon^{\mu\nu\rho\sigma} F_{\mu\nu} F_{\rho\sigma} = \frac{e^2}{2\pi^2}\, \mathbf{E}\cdot\mathbf{B}
```

微擾論上它來自三角圖；Fujikawa 的觀點更乾淨：路徑積分量測 𝒟ψ𝒟ψ̄ 在手徵轉動下有一個 Jacobian，用熱核 e^(−(iD̸)²/M²) 正規化後，正好得到 tr(γ₅ F F̃)。這同時是 Atiyah–Singer 指標定理：Dirac 算子的零模手徵數差，等於規範場的拓撲荷。實驗上的指紋是 π⁰ → γγ 的衰變率。

凝聚態與光子學的對應：Weyl 半金屬的負磁阻與手徵磁效應；光子晶體中的 Weyl 點（Lu 等人, Science 2015）。

### 習題（讀書筆記形式）

1. 用譜流論證，從頭推出 5.1 的結果，注意 Dirac 海的角色。
2. 從 Chern–Simons 作用量 (C/4π) ε A∂A 推出霍爾電流，驗證它恰好補償 5.2 的邊緣異常。
3. 在 1+1 維做 Fujikawa 計算：寫出手徵轉動下的 Jacobian，用熱核正規化求出異常項。
4. 寫一篇筆記：用 300 字向工程師解釋「為什麼第 4 講的疇壁一定有邊緣模」，同時用到第 4 講與本講的語言。

**閱讀順序**：David Tong《Gauge Theory》講義第 3 章〈Anomalies〉（3.1.1 譜流、3.2.2 Fujikawa 方法、3.3.1 指標定理；主教材，免費），anomaly inflow 見同講義 4.4.1 節；Callan & Harvey, Nucl. Phys. B 250, 427 (1985)；J. Harvey《TASI 2003 Lectures on Anomalies》（arXiv: hep-th/0509097）；Fujikawa & Suzuki《Path Integrals and Quantum Anomalies》（深入參考）；Nielsen & Ninomiya, Phys. Lett. B 130, 389 (1983)。

## 第 6 講：RCWA I——一維光柵的傅立葉模態法

核心結論：RCWA 在橫向用 Floquet 平面波展開（和 PWE 同一套數學），在縱向把每一層當成 z 方向均勻的波導，解一次本徵問題得到該層的所有模態，再用邊界條件把各層接起來。和 PWE 的差別是：本徵值不再是頻率，而是固定頻率下的縱向傳播常數。

**學習目標**：推導 TE 偏振下單層二元光柵的 RCWA；計算各繞射階的效率；用能量守恆與薄膜極限驗證程式。

### 6.1 幾何與 Floquet 展開

週期 Λ 的二元光柵（厚度 d，沿 x 週期、沿 y 均勻），上方入射區折射率 n\_I，下方基板 n\_II。平面波以入射角 θ 入射，TE 偏振（電場沿 y）。Floquet 定理給出第 m 繞射階的橫向波數：

```latex
k_{x,m} = k_0 \left( n_I \sin\theta - m \frac{\lambda}{\Lambda} \right), \qquad m = -M, \dots, M
```

共 N = 2M + 1 個諧波。光柵層內的場寫成 E\_y = Σ S\_m(z) exp(−i k\_(x,m) x)。

### 6.2 層內的本徵問題

以 z′ = k₀z 無因次化，K\_x = diag(k\_(x,m)/k₀)，E = ⟦ε⟧ 為介電常數的 Toeplitz 矩陣，波動方程化為二階常微分方程組：

```latex
\frac{d^2 \mathbf{S}}{dz'^2} = \mathbf{A}\,\mathbf{S}, \qquad \mathbf{A} = \mathbf{K}_x^2 - \mathbf{E} = \mathbf{W}\, \mathbf{Q}^2\, \mathbf{W}^{-1}
```

本徵向量 W 是該層的模態剖面，Q 的對角元（取實部為正、或虛部為正的分支）是模態的傳播常數。層內通解為 S(z′) = W \[e^(−Qz′) c⁺ + e^(Q(z′−d′)) c⁻\]。把 e^(Q(z′−d′)) 寫成以層底為參考點，是為了讓指數只會衰減、不會爆掉，這是第 8 講數值穩定性的伏筆。

對照 PWE：PWE 固定 k、解 ω；RCWA 固定 ω 與 k\_x、解 k\_z。兩者的「⟦ε⟧ 該怎麼取」是同一個問題。

### 6.3 邊界條件與繞射效率

入射區與基板是均勻介質，場就是 Rayleigh 展開。第 m 階的縱向波數為

```latex
k_{z,m}^{(\ell)} = \sqrt{(k_0 n_\ell)^2 - k_{x,m}^2}, \qquad \ell \in \{I, II\}
```

當平方根內為負，該階是衰減波（取虛部為正的分支）。在 z = 0 與 z = d 兩個介面匹配切向的 E\_y 與 H\_x，得到一組線性方程，解出反射振幅 R\_m 與穿透振幅 T\_m。TE 的繞射效率為

```latex
\mathrm{DE}_{r,m} = |R_m|^2\, \operatorname{Re}\!\left( \frac{k_{z,m}^{(I)}}{k_0 n_I \cos\theta} \right), \qquad \mathrm{DE}_{t,m} = |T_m|^2\, \operatorname{Re}\!\left( \frac{k_{z,m}^{(II)}}{k_0 n_I \cos\theta} \right)
```

**Rayleigh anomaly**：當某一階的 k\_z 恰好為零，該階正掠射離開、即將變成衰減波，效率曲線會出現尖銳的轉折。這就是 Wood anomaly 的一部分，用你的程式掃波長就能直接看到。

### 習題

1. 推導 6.2 的矩陣方程，並寫出 H\_x 用 W、Q 表示的形式。
2. 實作單層 TE 光柵。驗證（a）無損時 Σ DE\_r + Σ DE\_t = 1，誤差 < 10⁻¹⁰；（b）光柵退化為均勻薄膜時，與 Fresnel／薄膜干涉的解析式一致。
3. 掃描波長，畫出 0 階與 ±1 階效率；用 6.3 的條件預測 Rayleigh anomaly 的位置並與數值比對。
4. 畫出某一效率對 M 的收斂曲線。TE 應收斂很快，請記下來，第 7 講的 TM 會很不一樣。

**驗收標準**：能量守恆與薄膜極限兩項測試都納入自動化測試。

**閱讀**：Moharam、Grann、Pommet、Gaylord, JOSA A 12, 1068 (1995)（RCWA 的標準推導）；Moharam 等人, JOSA A 12, 1077 (1995)（增強穿透矩陣法）。

## 第 7 講：RCWA II——Li 的因式分解法則與二維週期結構

核心結論：TM 偏振與二維結構收斂慢，不是因為諧波不夠，而是把「兩個不連續函數的乘積」截斷錯了。Li 的三條定理告訴你每個乘積該用 Laurent 規則還是逆規則；用對了，金屬光柵的 TM 收斂可以從幾百階降到幾十階。

**學習目標**：證明並應用 Li 的因式分解定理；正確實作 TM 偏振；把 RCWA 推廣到二維交叉光柵，並理解任意形狀（如圓柱）的處理方法。

### 7.1 Li 的三條定理

設 h(x) = f(x) g(x)，三者都是週期函數，用截斷的傅立葉級數表示：

1. f 與 g 沒有**同時**的跳躍不連續 → Laurent 規則：⟦h⟧ = ⟦f⟧⟦g⟧。
2. f 與 g 有同時的跳躍，但乘積 h 連續 → 逆規則：⟦h⟧ = ⟦1/f⟧⁻¹⟦g⟧。
3. 有同時的跳躍且 h 也不連續 → 兩者都不正確，必須改寫問題。

應用到電磁場：在 x 方向的介面上，D\_x = εE\_x 連續、但 ε 與 E\_x 都跳，屬於定理 2；E\_z 連續而 ε 跳，D\_z = εE\_z 屬於定理 1。

### 7.2 TM 偏振的正確形式

對 1D 光柵的 TM（H\_y），層內本徵問題應寫成

```latex
\frac{d^2 \mathbf{U}}{dz'^2} = [\![1/\varepsilon]\!]^{-1} \left( \mathbf{K}_x\, [\![\varepsilon]\!]^{-1} \mathbf{K}_x - \mathbf{I} \right) \mathbf{U}
```

最前面那個 ⟦1/ε⟧⁻¹ 就是逆規則；1996 年以前的公式（包括第 6 講引用的 Moharam 1995）在這裡用 ⟦ε⟧，對介電光柵收斂慢，對金屬光柵幾乎不收斂。Lalanne & Morris、Granet & Guizal 在 1996 年從數值上發現這個修正，Li 同年給出數學證明。

### 7.3 二維交叉光柵

諧波變成二維指標 (m, n)，總數 N = (2M\_x + 1)(2M\_y + 1)，每層要解 2N × 2N 的本徵問題，成本是 O(N³) = O(M⁶)。M = 10 時 N = 441，本徵問題是 882 × 882；M = 20 時 N = 1681，已經是 3362 × 3362。這就是 2D RCWA 吃記憶體與時間的根源，也是 GPU 與演算法改良的空間。

介電張量的處理依形狀而定：

- **矩形柱**：Li（1997）的「先沿一個方向用逆規則、再沿另一個方向用 Laurent 規則」，對 x、y 兩個方向各組一次。
- **任意形狀（圓柱、斜邊）**：介面不再與座標軸對齊，要先在每一點定義介面法向量場，把 E 分解成法向（用逆規則）與切向（用 Laurent 規則）分量。這就是法向量場法（Schuster 等人, 2007）；S4 的 PolarizationDecomposition 選項即此類方法。

### 習題

1. 對兩個在同一點跳躍的步階函數，數值驗證定理 2：比較 ⟦f⟧⟦g⟧ 與 ⟦1/f⟧⁻¹⟦g⟧ 對 ⟦h⟧ 的截斷誤差隨 M 的變化。
2. 實作 1D TM，分別用舊公式與 7.2 的公式，對一個金屬（如 n = 0.2 + 3i 量級的複折射率）矩形光柵畫出 0 階反射率對 M 的收斂曲線。
3. 實作 2D 矩形柱光柵（Li 1997 規則），在單一偏振入射下驗證能量守恆。
4. 對圓柱陣列，比較「直接 FFT」與法向量場法的收斂速度。

**驗收標準**：1D TM 金屬光柵與文獻數值一致至 3 位有效數字；2D 矩形柱與 S4 的結果一致（第 8 講正式建立對照流程）。

**閱讀**：Li, JOSA A 13, 1870 (1996)；Lalanne & Morris, JOSA A 13, 779 (1996)；Granet & Guizal, JOSA A 13, 1019 (1996)；Li, JOSA A 14, 2758 (1997)（交叉光柵）；Schuster 等人, JOSA A 24, 2880 (2007)（法向量場）。

## 第 8 講：S 矩陣、數值穩定性與對 S4 的驗證

核心結論：用轉移矩陣（T 矩陣）串接多層時，衰減模的 e^(+Qz) 會指數爆炸，厚層或高階諧波必然失敗；改用散射矩陣（S 矩陣）只傳遞「入射 → 出射」的關係，所有指數都是衰減的，因此無條件穩定。求解器是否可信，則靠一套自動化的物理檢驗加上與 S4 的回歸比對來證明。

**學習目標**：推導 Redheffer 星積並實作多層 S 矩陣；辨識 RCWA 常見的數值陷阱；建立可重複的驗證流程。

### 8.1 為什麼 T 矩陣會失敗

T 矩陣把一個介面兩側的 (c⁺, c⁻) 線性連起來，串接時要同時處理 e^(−Qd) 與 e^(+Qd)。對高階的衰減模，Q 的實部很大，e^(+Qd) 在雙精度下溢位，或者讓矩陣的條件數大到結果全是雜訊。第 6 講把指數寫成以層底為參考點，就是避開正指數的第一步。

### 8.2 S 矩陣與 Redheffer 星積

S 矩陣把出射振幅寫成入射振幅的函數：

```latex
\begin{pmatrix} \mathbf{c}_1^- \\ \mathbf{c}_2^+ \end{pmatrix} = \begin{pmatrix} S_{11} & S_{12} \\ S_{21} & S_{22} \end{pmatrix} \begin{pmatrix} \mathbf{c}_1^+ \\ \mathbf{c}_2^- \end{pmatrix}
```

S₁₁、S₂₂ 是兩側的反射，S₂₁、S₁₂ 是兩個方向的穿透。兩個區塊 A（在前）與 B（在後）串接時，用 Redheffer 星積合併，各項就是多重反射的幾何級數求和：

```latex
\begin{aligned}
S_{11}^{AB} &= S_{11}^{A} + S_{12}^{A} \left( I - S_{11}^{B} S_{22}^{A} \right)^{-1} S_{11}^{B} S_{21}^{A} \\
S_{12}^{AB} &= S_{12}^{A} \left( I - S_{11}^{B} S_{22}^{A} \right)^{-1} S_{12}^{B} \\
S_{21}^{AB} &= S_{21}^{B} \left( I - S_{22}^{A} S_{11}^{B} \right)^{-1} S_{21}^{A} \\
S_{22}^{AB} &= S_{22}^{B} + S_{21}^{B} \left( I - S_{22}^{A} S_{11}^{B} \right)^{-1} S_{22}^{A} S_{12}^{B}
\end{aligned}
```

每一層的 S 矩陣只含 e^(−Qd)，所以整個流程沒有正指數。實作時用 solve 取代顯式求反；Rumpf（2011）建議把每層夾在零厚度的「間隙介質」之間，讓各層的 S 矩陣對稱、程式更簡潔。

### 8.3 數值陷阱清單

| 陷阱 | 症狀 | 處理 |
| --- | --- | --- |
| 平方根分支選錯 | 衰減模被當成增長模，結果爆掉 | 統一取 Im(k\_z) ≥ 0，傳播模取 Re(k\_z) > 0 |
| 均勻層的簡併本徵值 | 本徵向量任意混合、W 病態 | 均勻層直接用解析平面波，不呼叫本徵解法 |
| Rayleigh anomaly 附近 k\_z → 0 | 矩陣近奇異 | 輕微偏移波長或加極小損耗做檢查 |
| 截斷方式不一致 | 與參考程式在同一 M 下對不上 | 比較收斂後的值，而非同一 M 的值 |
| 顯式矩陣求反 | 精度損失 | 一律用 LU 分解求解 |

### 8.4 驗證流程

每次提交都自動跑以下測試：

1. **能量守恆**：無損結構的總效率等於 1。
2. **互易性**：交換入射與出射方向，對應的穿透係數相等。
3. **對稱性**：對稱結構在 θ 與 −θ 下得到鏡像的繞射效率。
4. **極限**：均勻薄膜對解析解；單層 TE 對第 6 講結果。
5. **收斂**：關鍵輸出隨 M 的變化曲線，存成圖檔隨版本追蹤。
6. **回歸**：與 S4 在同一組結構上比對收斂值。

Fan & Joannopoulos（2002）分析的光子晶體平板是很好的基準（用 S4 計算同一結構作對照）：介電平板上打方晶格空氣孔，正入射穿透譜會出現尖銳的 Fano 型導模共振。共振的位置與線寬對截斷與因式分解都很敏感，是比平滑光譜嚴格得多的考驗。

### 習題

1. 由多重反射的幾何級數推導 Redheffer 星積的四個式子。
2. 用 T 矩陣與 S 矩陣分別計算一個厚度遞增的光柵，找出 T 矩陣開始失敗的厚度與諧波數。
3. 把 8.4 的六項測試寫成 pytest，納入持續整合（GitHub Actions）。
4. 重現光子晶體平板的穿透譜，與 S4 比對共振頻率與 Q 值。

**驗收標準**：六項測試全數通過；光子晶體平板共振頻率與 S4 相差 < 10⁻³（相對）。

**閱讀**：Li, JOSA A 13, 1024 (1996)（S 矩陣與 R 矩陣演算法的比較）；Rumpf, PIER B 35, 241 (2011)；Liu & Fan, Comput. Phys. Commun. 183, 2233 (2012)（S4）；Fan & Joannopoulos, PRB 65, 235112 (2002)。

## 第 9 講：求解器軟體工程——C++ 核心、pybind11、測試與發佈

核心結論：Python 原型是「正確答案的神諭」，C++ 核心只負責更快地給出同一個答案；兩者由同一套 pytest 驅動，誤差要求到 10⁻¹²。EDA 公司看重的不只是你算得對，而是你的程式能被別人安裝、測試、維護。

**學習目標**：設計分層的求解器架構；用 pybind11 做零拷貝的 NumPy 介面；建立 CMake 建置、持續整合與 wheel 發佈流程。

### 9.1 架構

1. **幾何層**：週期、層的厚度與材料、圖形（矩形、圓、多邊形）。
2. **傅立葉係數引擎**：解析式或 FFT，加上 Laurent／逆規則／法向量場的選擇。
3. **層求解器**：組本徵矩陣、求 W 與 Q；均勻層走解析捷徑。
4. **S 矩陣組裝**：逐層計算、星積串接。
5. **輸出**：繞射效率、場重建、（第 11 講）梯度。

每一層都有獨立的單元測試，並可個別替換。例如之後把第 3 層換成 GPU 版本，其他層不動。

### 9.2 技術選擇

- **線性代數**：用 Eigen（純標頭檔、易於建置）起步；本徵問題是效能瓶頸，之後可提供 LAPACK（zgeev）後端。不建議自己寫本徵解法，那是另一個專業領域，也不是 EDA 公司要考你的東西。
- **綁定**：pybind11，搭配 Eigen 型別轉換與 NumPy 緩衝區協定，避免大矩陣複製。
- **建置**：CMake + scikit-build-core；用 cibuildwheel 在 GitHub Actions 上為 Linux、Windows、macOS 產生 wheel。
- **專案位置**：repo 放在 OneDrive 等同步資料夾之外，避免建置時檔案被鎖。

### 9.3 工作流程

- Python 參考實作與 C++ 核心並存；同一組 pytest 參數化地跑兩者，差異 > 10⁻¹² 即失敗。
- 每個功能用 PR 合併，即使只有你一個人：PR 描述就是你的設計紀錄，也是面試時可以展示的東西。
- 文件用 MkDocs 或 Sphinx；範例庫（光柵、光子晶體平板、Chern 數）同時是部落格文章的素材。
- 授權選 MIT 或 BSD；README 明確寫出「僅依據公開文獻實作」。

### 9.4 效能剖析

在最佳化之前先量測。2D RCWA 的時間通常絕大部分花在每層的本徵分解，其次是星積中的線性求解；傅立葉係數與場重建相對便宜。用 perf 或 VTune 看熱點，並把「每層時間 vs 諧波數」畫成 log–log 圖，確認是否符合 O(N³)。

### 習題

1. 把第 6 講的 1D TE 求解器移植到 C++，透過 pybind11 呼叫，與 Python 版本逐位比對。
2. 量測 NumPy 陣列傳入 C++ 時是否發生複製，並改成零拷貝。
3. 建立 GitHub Actions：每次推送都跑測試，標籤推送時建置三平台 wheel 並上傳到 TestPyPI。
4. 對 2D RCWA 做效能剖析，畫出時間對諧波數的關係，並與 Python 版比較加速倍率。

**驗收標準**：在一台全新的機器上 pip install 後，範例庫全部能直接執行。

**閱讀**：pybind11 官方文件（Eigen 與 NumPy 章節）；scikit-build-core 文件；Stanford CS106L（現代 C++ 的慣用寫法）。

## 第 10 講：FDTD 與 GPU 加速

核心結論：FDTD 是在 Yee 交錯網格上做蛙跳時間積分的模板（stencil）計算，每個網格點的運算量極少、記憶體存取極多，所以在 GPU 上是記憶體頻寬受限的問題。加速的上限由頻寬決定，衡量成果要看「達到峰值頻寬的幾成」，而不只是比 CPU 快幾倍。

本講對應你在公司的 FDTD 工作。公開作品只放一般性的玩具實作與不含機密的指標，公司的程式碼與參數不出公司。

**學習目標**：推導 Yee 格式、CFL 條件與數值色散；實作 CPML 吸收邊界；用 roofline 模型分析並優化 GPU 核心；用 FDTD 與 PWE 交叉驗證。

### 10.1 Yee 格式與穩定性

E 與 H 在空間上錯開半格、時間上錯開半步。以 1D 為例：

```latex
H_y^{n+1/2}(i+\tfrac12) = H_y^{n-1/2}(i+\tfrac12) - \frac{\Delta t}{\mu \Delta x}\left[ E_z^{n}(i+1) - E_z^{n}(i) \right], \qquad E_z^{n+1}(i) = E_z^{n}(i) - \frac{\Delta t}{\varepsilon \Delta x}\left[ H_y^{n+1/2}(i+\tfrac12) - H_y^{n+1/2}(i-\tfrac12) \right]
```

（符號依傳播方向與場分量的取向而定，推導時請從旋度方程自行確認。）三維的 CFL 穩定條件為

```latex
\Delta t \le \frac{1}{c\,\sqrt{\Delta x^{-2} + \Delta y^{-2} + \Delta z^{-2}}}
```

數值色散讓高頻成分傳得太慢。常見的經驗法則是在最高折射率材料中每個波長取 10–20 格以上，但最終要以收斂測試為準。對 CIS 像素這類亞波長結構，真正決定網格的往往是幾何細節而非波長。

### 10.2 邊界、光源與材料

- **吸收邊界**：CPML（Roden & Gedney, 2000），以拉伸座標加遞迴卷積實作複頻移 PML，與宿主介質無關，對衰減波也有良好吸收。
- **週期邊界**：像素陣列用週期或 Bloch 邊界；斜入射的寬頻模擬需要特別處理（例如分裂場法），是 FDTD 相對於 RCWA 的弱點。
- **光源與量測**：總場／散射場（TF/SF）注入平面波；用寬頻脈衝加上執行中的 DFT 監測器，一次模擬得到整段光譜。
- **色散材料**：金屬與半導體用 Drude／Lorentz 模型，以輔助微分方程（ADE）加入更新式。

### 10.3 GPU 優化：roofline 思維

每個網格點的更新只讀寫少數幾個浮點數、做十幾次運算，算術強度很低，所以效能上限 ≈ 記憶體頻寬 ÷ 每格位元組數。優化的順序：

1. **資料佈局**：SoA（每個場分量一個陣列），讓相鄰執行緒讀相鄰位址，達成合併存取。
2. **核心融合**：減少全域記憶體的往返次數；把係數陣列壓縮（例如用材料索引查表取代每格存係數）。
3. **共享記憶體分塊**：對 stencil 的鄰點重複讀取有幫助，但在現代 GPU 的 L1／L2 快取下收益常常有限，要量測才知道。
4. **精度**：場更新多數情況 float32 足夠；長時間累加的 DFT 監測器要評估是否需要 float64。

報告指標建議用每秒更新的網格數（Mcells/s）與「實測頻寬 ÷ 峰值頻寬」，再用 Nsight Compute 佐證。

### 習題

1. 寫 1D FDTD 加 CPML，量測邊界反射係數對 PML 厚度的關係。
2. 量測數值色散：比較模擬的相速度與理論值，對每波長格數作圖。
3. 寫 2D TM FDTD，以 Bloch 週期邊界與隨機點源激發光子晶體，用 harminv 擷取共振頻率，與第 2 講的 PWE 能帶交叉驗證。
4. 把習題 3 的更新核心移植到 GPU（CUDA 或 CuPy），畫出 roofline，報告達到峰值頻寬的比例。

**驗收標準**：FDTD 擷取的能帶頻率與 PWE 相差 < 1%；GPU 核心達到峰值記憶體頻寬的 60% 以上。

**閱讀**：Taflove & Hagness《Computational Electrodynamics: The Finite-Difference Time-Domain Method》（Yee 格式、數值穩定性與 PML 相關章節）；Roden & Gedney, Microwave Opt. Technol. Lett. 27, 334 (2000)；Oskooi 等人, Comput. Phys. Commun. 181, 687 (2010)（Meep）；NVIDIA《CUDA C++ Best Practices Guide》。

## 第 11 講：可微分求解器與 adjoint 逆向設計

核心結論：adjoint 方法只需多解一次「伴隨問題」，就能得到目標函數對任意多個設計參數的梯度，成本與參數數量無關；這是梯度式逆向設計能處理上萬個自由度、而遺傳演算法不行的原因。對 RCWA 而言，難點在於如何穩定地對本徵分解求導。

**學習目標**：推導線性系統的 adjoint 梯度；實作本徵分解的導數並處理簡併；用梯度檢查驗證；完成一個可公開的逆向設計示範。

### 11.1 Adjoint 的一般形式

設正向問題為 A(p) x = b(p)，目標為實數 f(x, p)。引入伴隨變數 λ 滿足 Aᵀλ = (∂f/∂x)ᵀ，則

```latex
\frac{df}{dp_j} = \frac{\partial f}{\partial p_j} - \boldsymbol{\lambda}^{\mathsf{T}} \left( \frac{\partial A}{\partial p_j}\, \mathbf{x} - \frac{\partial \mathbf{b}}{\partial p_j} \right)
```

正向一次、伴隨一次，之後對每個 p\_j 只需一個便宜的矩陣–向量乘積。場是複數、目標是實數時，要用 Wirtinger 導數並取實部，這是最常出錯的地方，務必用梯度檢查把關。

### 11.2 對 RCWA 求導：本徵分解的導數

RCWA 的每一層經過 A = W Λ W⁻¹。對非簡併的本徵值，一階微擾給出

```latex
d\Lambda_{ii} = \left( W^{-1}\, dA\, W \right)_{ii}, \qquad dW = W \left( F \circ \left( W^{-1}\, dA\, W \right) \right), \quad F_{ij} = \begin{cases} \dfrac{1}{\lambda_j - \lambda_i} & i \ne j \\ 0 & i = j \end{cases}
```

當 λ\_i ≈ λ\_j（均勻層、對稱結構都會出現），F 發散。實務做法是 Lorentzian 展寬 F\_ij ≈ (λ\_j − λ\_i)/((λ\_j − λ\_i)² + δ²)，或把簡併子空間當成整體處理。開源的 grcwa 與 TORCWA 分別示範了自動微分與這類處理，值得逐行對照。

兩條實作路線：

- **自動微分**：在 JAX 或 PyTorch 上重寫 Python 版 RCWA，讓框架追蹤導數；開發快，適合研究與原型。
- **手寫 adjoint**：在 C++ 核心內直接推導並實作梯度；速度與記憶體最佳，是商用工具的做法，也最能展現你懂底層。

### 11.3 拓撲最佳化的規則化

像素化的自由形狀設計，需要三個步驟把連續密度變成可製造的結構：密度濾波（控制最小特徵尺寸）、tanh 投影（逐步推向 0／1）、二值化後的最終驗證。投影的陡度要分階段提高，否則最佳化會卡在灰階解。

### 11.4 與遺傳演算法的比較

遺傳演算法不需要梯度、能跳出局部極值，但評估次數隨參數數量快速增加，適合十幾個參數的形狀最佳化。adjoint 適合成千上萬個像素的自由形狀，但只保證局部最優，通常需要多個隨機初始點。實務上常先用 adjoint 找出好的拓撲，再用少量參數做穩健性微調。

### 習題

1. 推導 11.1，並寫出場為複數、目標為 |場|² 時的完整形式。
2. 用 JAX 實作 1D RCWA（TE），對光柵寬度與厚度求梯度；用中央差分做梯度檢查，相對誤差 < 10⁻⁶。
3. 實作 11.2 的本徵導數，在均勻層與對稱光柵上重現簡併造成的發散，並用 Lorentzian 展寬修正。
4. 示範專案：最佳化一個週期性超表面光束偏折器（metagrating），最大化目標波長的 +1 階繞射效率；最終結構用獨立求解器（S4 或 FDTD）驗證。
5. 對同一問題用遺傳演算法、限制在 10 個參數以內，比較收斂所需的模擬次數與最終效率。

**驗收標準**：梯度檢查通過；最佳化結果在獨立求解器中重現，效率差異 < 1%。

**閱讀**：Molesky 等人, Nature Photonics 12, 659 (2018)（逆向設計綜論）；Hughes 等人, ACS Photonics 5, 4781 (2018)（adjoint 推導，含非線性推廣）；Jin 等人, ACS Photonics 7, 2350 (2020)（grcwa，以 autograd 微分的 RCWA）；Sell 等人, Nano Letters 17, 3752 (2017)（自由形狀 metagrating）；Kim & Lee, Comput. Phys. Commun. 282, 108552 (2023)（TORCWA）；Hammond 等人, Optics Express 30, 4467 (2022)（Meep 的 adjoint 與製程約束）。

## 附錄：里程碑檢核表與參考文獻

每完成一項就打勾；一季內沒有任何一項可以打勾，代表範圍需要縮小。

**第一階段（0–3 個月）**

- [ ] 確認勞動契約中的兼職、智財與競業條款
- [ ] 建立公開 repo（放在同步資料夾之外），設定 pytest 與 GitHub Actions
- [ ] 1D PWE 對 Rytov 解析解（第 2 講習題 1）
- [ ] 2D TM 能帶對 MPB，誤差 < 0.5%
- [ ] QWZ 模型的 Chern 數相圖
- [ ] 旋磁光子晶體能隙以下 Chern 數總和的絕對值為 1
- [ ] 疇壁超晶胞的邊緣模
- [ ] 技術筆記 ×2

**第二階段（3–9 個月）**

- [ ] 1D TE RCWA，通過能量守恆與薄膜極限測試
- [ ] 1D TM 金屬光柵，Li 法則前後的收斂比較
- [ ] 多層 S 矩陣與六項驗證測試
- [ ] 光子晶體平板穿透譜對 S4
- [ ] C++ 核心 + pybind11，pip 可安裝
- [ ] 公司 FDTD GPU 工作的可公開指標（Mcells/s、頻寬效率）

**第三、四階段（9–24 個月）**

- [ ] JAX 版可微分 RCWA，梯度檢查通過
- [ ] 本徵導數的簡併處理
- [ ] metagrating 逆向設計示範，經獨立求解器驗證
- [ ] 研討會口頭或海報發表
- [ ] 以作品集應徵 EDA／光學模擬職缺

### 參考文獻

以下書目與卷頁已於 2026-10-29 對照出版社、arXiv、作者或機構頁面逐一查核，主要查核頁面列於文末「查核來源」。

**教科書與講義**

- Joannopoulos、Johnson、Winn、Meade《Photonic Crystals: Molding the Flow of Light》第 2 版（2008）
- Sakoda《Optical Properties of Photonic Crystals》
- Taflove & Hagness《Computational Electrodynamics: The FDTD Method》
- Asbóth、Oroszlány、Pályi《A Short Course on Topological Insulators》（arXiv:1509.02295）
- David Tong《Gauge Theory》講義，劍橋大學
- Fujikawa & Suzuki《Path Integrals and Quantum Anomalies》

**PWE 與因式分解**：Ho、Chan、Soukoulis（1990）；Li（1996, 1997）；Lalanne & Morris（1996）；Granet & Guizal（1996）；Johnson & Joannopoulos（2001，MPB）；Schuster 等人（2007）；Sözüer、Haus、Inguva（1992）；Farjadpour 等人（2006）。

**拓撲光子學**：Fukui、Hatsugai、Suzuki（2005）；Qi、Wu、Zhang（2006）；Xiao、Chang、Niu（2010）；Haldane & Raghu（2008）；Wang 等人（2008, 2009）；Lu、Joannopoulos、Soljačić（2014）；Lu 等人（2015，光子 Weyl 點）。

**量子異常**：Nielsen & Ninomiya（1983）；Callan & Harvey（1985）；Harvey, hep-th/0509097。

**RCWA 與 S 矩陣**：Moharam 等人（1995 兩篇）；Li（1996，S 矩陣）；Rumpf（2011）；Liu & Fan（2012，S4）；Fan & Joannopoulos（2002）。

**FDTD 與 GPU**：Roden & Gedney（2000）；Oskooi 等人（2010，Meep）；NVIDIA CUDA C++ Best Practices Guide。

**逆向設計**：Molesky 等人（2018）；Hughes 等人（2018）；Sell 等人（2017）；Kim & Lee（2023，TORCWA）；Jin 等人（2020，grcwa）；Hammond 等人（2022）。

### 查核來源

- [Wang 等人 2008（arXiv:0712.1776，全文）](https://ar5iv.arxiv.org/html/0712.1776)：YIG 參數、能隙、介面結構
- [Ho、Chan、Soukoulis 1990（全文 PDF）](http://repository.ust.hk/ir/bitstream/1783.1-43967/1/PhysRevLett.65.3152.pdf)：確認該文主題為鑽石晶格 3D 能隙
- [Tong《Gauge Theory》第 3 章](https://www.damtp.cam.ac.uk/user/tong/gaugetheory/3anom.pdf)：1+1 與 3+1 維異常公式、章節編號
- [MPB Python 教學](https://mpb.readthedocs.io/en/latest/Python_Tutorial/)：方晶格 ε = 12、R = 0.2a 的 TM 能隙
- [Fukui、Hatsugai、Suzuki（arXiv:cond-mat/0503172）](https://arxiv.org/abs/cond-mat/0503172)
- [Xiao、Chang、Niu（arXiv:0907.2021）](https://arxiv.org/abs/0907.2021)
- [Asbóth 等人（arXiv:1509.02295）](https://arxiv.org/abs/1509.02295)
- [Haldane & Raghu（arXiv:cond-mat/0503588）](https://arxiv.org/abs/cond-mat/0503588)
- [Lu 等人 2015，Weyl 點（MIT DSpace）](https://dspace.mit.edu/handle/1721.1/98069)
- [Harvey, TASI 2003（arXiv:hep-th/0509097）](https://arxiv.org/abs/hep-th/0509097)
- [Callan & Harvey 1985（Princeton）](https://collaborate.princeton.edu/en/publications/anomalies-and-fermion-zero-modes-on-strings-and-domain-walls/)
- [Moharam 等人 1995（UCF STARS）](https://stars.library.ucf.edu/scopus1990/1846)
- [Johnson & Joannopoulos 2001，MPB（全文）](https://math.mit.edu/%7Estevenj/papers/JohnsonJo01.pdf)
- [Farjadpour 等人 2006（全文）](https://math.mit.edu/%7Estevenj/papers/FarjadpourRo06.pdf)
- [Oskooi 等人 2010，Meep（全文）](https://math.mit.edu/~stevenj/papers/OskooiRo10.pdf)
- [Meep 致謝頁（Hammond 2022 引用格式）](https://meep.readthedocs.io/en/latest/Acknowledgements/)
- [Rumpf 2011 講義](https://empossible.net/wp-content/uploads/2020/05/Lecture-Scattering-Matrices-for-Semi-Analytical-Methods.pdf)
- [Fan & Joannopoulos 2002（全文）](https://uio.no/studier/emner/matnat/its/TEK5400/h18/TEK%209400%202018/papers/fan_prb_65_235112_2002.pdf)
- [Molesky 等人 2018（Princeton）](https://collaborate.princeton.edu/en/publications/inverse-design-in-nanophotonics/)
- [Hughes 等人 2018（arXiv:1811.01255）](https://arxiv.org/abs/1811.01255v1)
- [TORCWA（GitHub）](https://github.com/kch3782/torcwa)
- [grcwa（PyPI）](https://pypi.org/project/grcwa)
- [Fujikawa & Suzuki（Princeton 圖書館目錄）](https://catalog.princeton.edu/catalog/SCSB-8743797)
