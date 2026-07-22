# RCC 免疫治疗中突变负荷与 HLA-I 纯合性的治疗背景依赖性：多队列再分析及等位基因层面肽审计

## 摘要

### 背景

突变负荷与人类白细胞抗原 I 类（HLA-I）多样性分别描述肿瘤抗原性的不同环节，但二者与肾细胞癌（RCC）免疫检查点抑制剂（ICI）疗效的关系并不一致。本研究考察 lenvatinib 联合 pembrolizumab 队列中观察到的突变-HLA 现象是否可被统计学识别、能否在其他接受治疗的 RCC 队列中复现，以及是否得到等位基因层面肽预测的支持。

### 方法

基于 Lee 等报道的 24 例 RCC 队列，使用原文公开的 MAF 与 HLA 补充表及经独立核对的无进展生存期（PFS）终点重建分析。突变负荷以蛋白改变突变数表示，不将其称为临床 TMB。发现队列采用 Kaplan-Meier、log-rank、12 个月限制平均生存时间（RMST）、固定组大小置换、连续 Cox 与阈值敏感性分析。随后在 Miao/DFCI 的 35 例 nivolumab 治疗 ccRCC 队列中开展 mutation-only 分析，并对 Braun、JAVELIN Renal 101、Chowell、Nature Cancer 和 MSK 数据进行具有直接性分层的外部审计。MSK 分析将 ICI 结局与 MSK50K HLA 字段按样本 ID 精确匹配，对 20 个条件检验和 10 个交互检验分别校正多重比较。最后，用 MHCflurry 2.2.1 对全部 63 个归档保留的、位置一致的 HLA-A*02:01 与 partner allele 9-mer 事件进行重新评分。

### 结果

发现队列中 3 例 HLA-A 纯合患者均高于官方突变数中位数且均发生 PFS 事件。与其余 21 例相比，其 PFS 差异的 log-rank P=0.00045，固定组大小置换 P=0.0094，12 个月 RMST 差异为 -4.40 个月。但“高突变数/HLA 纯合”状态与 HLA-A 纯合完全重合，低突变数/纯合单元为空；连续突变数与 PFS 无关（每 IQR HR=1.48，95%CI 0.72-3.05，P=0.283）。在 Miao/DFCI 队列中，较高蛋白改变突变数与较短 PFS 相关（每 IQR HR=2.38，95%CI 1.12-5.04，P=0.024），但较大的 nivolumab 和 avelumab-axitinib 队列未复现简单的突变-HLA 关联。MSK 高 TMB ccRCC 中任一 HLA-I 位点纯合与较短 OS 相关（HR=4.83，95%CI 1.89-12.33，FDR=0.0196），但最强的 TMB×HLA 交互项未通过校正（FDR=0.293）。在完整肽事件宇宙中，partner allele 在 29/63 个事件中同时优于 HLA-A*02:01（双侧二项检验 P=0.615）；分配模式主要与肽 C 端残基类别相关。

### 结论

RCC 中突变负荷和 HLA-I 纯合性与结局的关系具有队列和治疗背景依赖性。发现队列复合状态不能被识别为突变-HLA 交互作用，MSK 结果是后验条件关联而非直接复现。等位基因层面的预测支持 C 端残基依赖的互补性，而不支持第二个 HLA-A 等位基因的全局“救援”。本研究界定了后续前瞻性免疫基因组研究需验证的边界，不能据此建立临床生物标志物。

**关键词：** 肾细胞癌；免疫检查点抑制剂；酪氨酸激酶抑制剂；HLA-I；突变负荷；HLA 纯合性；新抗原预测

## 引言

晚期 ccRCC 的一线治疗已广泛采用 ICI 联合 TKI 方案[2,3,22,23]。尽管总体疗效显著，不同患者获得的获益深度和持续时间差异很大。与高突变实体瘤不同，PD-L1、总突变负荷和单一免疫特征在 RCC 中均未形成稳定、可迁移的预测指标[4-8]。

突变负荷与 HLA-I 基因型不是同一层级的变量。突变提供候选抗原底物，而 HLA 基因型限制哪些肽能够被呈递。泛癌研究提示 HLA-I 全杂合性及较高 HLA evolutionary divergence（HED）与 ICI 结局相关[11,12]；但 ccRCC 具有较低突变负荷、显著的拷贝数与微环境异质性，以及受到治疗方案影响的免疫反应，因此不能默认上述关系在 RCC 中方向一致[4,18,19]。

本研究来源于 Lee 等的 24 例 lenvatinib-pembrolizumab 队列。原文已经报道 mean HED 与 PFS 和缓解持续时间相关，也报告 4 例任一 HLA-I 位点纯合患者呈不良趋势[1]。因此，本研究的贡献不是再次“发现 HLA 多样性”。需要回答的是更窄的问题：该队列中的位点特异性突变-HLA 模式能否被统计学区分；它在其他接受治疗的 RCC 队列中是否具有一致方向；以及配对肽预测是否支持一般性的呈递救援机制。

我们据此进行来源分层再分析：发现队列使用原文官方 MAF 和 HLA 表重建；外部证据同时纳入支持和不支持的队列；旧的 29 个方向筛选肽事件被完整的 63 个配对事件宇宙替代。本文定位为多队列的边界审计，而非验证性生物标志物研究。

## 方法

### 研究设计与证据层级

本研究为公开去标识化数据的二次分析。证据按以下层级处理：原始论文官方补充表与链接结局；可独立访问的接受治疗队列患者级数据；定义不完全一致的代理分析；以及不含患者结局的生物学背景。发现关联、外部关联、交互检验、肽预测和免疫肽组背景不互相替代。

### 发现队列重建

发现队列为 Lee 等报道的 24 例 lenvatinib-pembrolizumab 治疗 RCC 患者（NCT02501096；PRJNA610643）[1]。体细胞变异来自 Supplementary Data 1 的 `MAF exome sequencing` 表，HLA-I 基因型和发表的 mean HED 来自 Supplementary Data 2。按患者、基因组坐标、等位基因和变异分类去重后，官方表共有 1,296 条蛋白改变记录。每例蛋白改变突变数作为突变负荷指标；由于官方补充表不提供可调用区域分母，该指标不称为 TMB。

PFS 时间和删失状态经原文 Kaplan-Meier risk table 与 swimmer plot 核对后锁定为 16 个事件、8 个删失。历史 `Event` 字段与 10 例记录冲突，未用于本研究。HLA-A、-B、-C 纯合定义为原文发表的同一位点两条等位基因一致。官方与内部 HLA 基因型在 21/24 例三个位点均一致，主要分析均使用官方基因型。

### 统计分析

主要描述状态定义为突变数不低于队列中位数且任一 HLA-I 位点纯合。在该 24 例队列中，该状态恰好与 HLA-A 纯合完全相同。PFS 使用 Kaplan-Meier 与双侧 log-rank 比较，并计算 12 个月 RMST 差异。固定组大小置换将观察到的 log-rank 统计量与 24 例中所有 3 例组合比较；由于该状态是在观察数据后定义，置换并不能校正完整的模型选择过程。

蛋白改变突变数进一步按 IQR 作为连续变量进行单变量 Cox 分析，并使用 scaled Schoenfeld residuals 检查比例风险。敏感性分析以突变数上三分位数替代中位数，并单独考察任一位点纯合。由于没有 HLA-A 纯合患者位于中位数以下，发现队列不拟合多变量交互模型。

### 外部队列

Miao/DFCI 公共队列 `ccrcc_dfci_2019` 中，35 例 nivolumab 治疗 ccRCC 患者有 PFS、RECIST 和 2,824 条蛋白改变突变记录[5]。每例蛋白改变突变数以秩变量表示并按 IQR 缩放；进一步调整一线治疗状态和 log2 nivolumab 剂量，采用二次项检验非线性并进行逐例删除诊断。35 例拟合样本的三个模型变量均无缺失，定义见补充表 S7。该数据不提供 HLA、HED、HLA LOH 或个体化呈递数据，因此只能检验 mutation-only 关联。

外部 RCC 队列不进行合并 meta-analysis，因为治疗、终点、测序平台与 HLA 构建不同。Braun nivolumab 数据包含 261 例患者的 PFS、OS、突变数、新抗原数和 HLA-I 纯合性[4]。JAVELIN Renal 101 源数据包含 333 例 avelumab-axitinib 治疗患者的每 Mb 非同义变异、HLA-A subtype 数和 PFS[7]。Chowell 多 ICI RCC 表含 91 例患者的 TMB、HED、HLA LOH、PFS 和 OS[14]。Nature Cancer 源数据提供 128 例 IO/IO 患者按反应类别分层的 mean HED[8]。SNiP-RCC 仅作为已发表的 HLA-only 背景结果引用[9]。

MSK 分析将 `tmb_mskcc_2018` ICI 结局与 `msk_impact_50k_2026` HLA 字段按样本 ID 精确匹配[10,15]，每例仅保留一个样本，得到 138 例 RCC、51 个 OS 事件和 115 例 ccRCC。考察 HLA-A 纯合、任一 HLA-I 纯合、HLA-A LOH、任一 HLA-I LOH 和 HED-A=0。对 2 个 TMB 阈值、2 个组织学范围和 5 个 HLA 定义进行 20 个条件检验，并在该家族内进行 BH 校正；5 个连续与 5 个二分类 TMB×HLA 交互项构成独立的 10 项校正家族。模型调整年龄、性别及联合 ICI 与单药治疗。MSK50K 采用 CC BY-NC-ND 4.0 许可，因此不随稿重新发布修改后的患者级表，仅提供研究 ID、代码和聚合输出。

### HLA-A*02:01 配对肽审计与免疫肽组背景

归档肽表中有来自 7 例 HLA-A*02:01 杂合患者的 63 个归档保留、位置一致、突变来源 C 端 9-mer 事件，其中 34 个原始标记为 HLA-A*02:01，29 个标记为 partner allele。29 个 partner 事件不能单独推断，因为它们已按待检验方向选择。全部 63 个事件使用 MHCflurry 2.2.1 和 presentation model release 2.2.0 以原始 N/C flank 对两条 HLA-A 等位基因重新评分[16]；LEN07 使用发表的 HLA-A*34:01。partner win 需同时满足更低 affinity、更低 affinity percentile、更高 presentation score 和更低 presentation percentile。总 partner win 比例以双侧二项检验比较 0.5；C 端分为 basic（K/R/H）、aromatic（F/W/Y）、aliphatic hydrophobic（L/I/V/M/A）和其他，采用 Fisher 精确检验并在 4 类内进行 BH 校正。事件嵌套于 7 例患者，以上 P 值仅为事件层面的描述性检验，不能作为患者层面推断。

PXD017149 的 1,457 条全部 ccRCC HLA-I 配体和 504 条肿瘤特异性配体用于计算末端序列和 C 端残基类别构成[17]。该资源不含本研究患者突变、HLA 配对或结局，仅作为真实配体背景。PASS Mutect2 SNV 的 COSMIC SBS v3.2 严格 bootstrap 审计显示 3 个历史核心病例的 SBS2 和 SBS2+SBS13 检出率均为 0/100，因此 APOBEC 标签不纳入主分析。

## 结果

### 官方来源重建后，发现队列的暴露定义发生改变

官方 Lee MAF 的每例蛋白改变突变数中位数为 48.5（IQR 27.2），与原文报道一致[1]。该指标与历史服务器保留记录数不可互换，两者患者排序相关性较弱（Spearman rho=0.166，P=0.438）。因此，后续仅使用官方突变数与发表 HLA 基因型配对（Figure 1）。3 例 HLA-A 纯合患者为 EIS01、LEN08 和 LEN15；另有 1 例在 HLA-B 与 HLA-C 纯合。

![Figure 1](figures/Figure1.png)

### 发现状态对应早进展，但不能识别为突变-HLA 交互作用

3 例 HLA-A 纯合患者均高于官方突变数中位数，且均发生 PFS 事件。其与其余 21 例的 PFS 差异为 log-rank P=0.00045；固定 3 例分组置换 P=0.0094；12 个月 RMST 差异为 -4.40 个月（Figure 2）。这些结果只说明这 3 条观察到的 PFS 记录相对于其他固定大小的三人组较早，并未校正 HLA 位点、阈值及复合状态的后验选择。

更重要的是，该复合状态与 HLA-A 纯合完全相同，低突变数/HLA-A 纯合单元为空。连续突变数与 PFS 无显著关联（每 27.2 个突变 HR=1.48，95%CI 0.72-3.05，P=0.283）；采用上三分位阈值时核心组由 3 例降至 2 例。任一 HLA-I 位点纯合单独定义时包含 4 例，PFS 差异不显著（log-rank P=0.253）。因此，该结果只能作为位点特异病例观察，不能拆分为独立的突变和 HLA 效应。

![Figure 2](figures/Figure2.png)

### Miao/DFCI 支持 mutation-only 关联，但外部 RCC 队列总体不一致

Miao/DFCI 中突变秩 IQR 为 0.486，突变数 IQR 为 34。较高突变秩与较短 PFS 相关（每 IQR HR=2.38，95%CI 1.12-5.04，P=0.024；比例风险检验 P=0.834）。原始突变数模型方向相近（每 34 个突变 HR=2.43，95%CI 1.32-4.48，P=0.004）；调整一线治疗和剂量后仍存在关联（HR=2.47，P=0.017）。35 次逐例删除均保持 HR>1（2.04-2.78），34/35 次名义 P<0.05。中位数分组的 log-rank P=0.133，故其仅用于描述（Figure 3）。

![Figure 3](figures/Figure3.png)

该方向并未在所有外部队列中重复。Braun nivolumab ccRCC 中，高突变负荷、高新抗原负荷、HLA-I 纯合及其组合均与 PFS/OS 无关；高突变负荷加任一位点纯合的 PFS P=0.651、OS P=0.608。JAVELIN avelumab-axitinib 高负荷亚组中，HLA-A 单一 subtype 与多个 subtype 的 PFS 无差异（调整 HR=0.91，95%CI 0.45-1.81，P=0.780）。Chowell RCC 子集的 high-TMB/low-HED 与 PFS、OS 均无关；Nature Cancer IO/IO 中 mean HED 在反应类别间无差异。相反，SNiP-RCC 报道较高 HED-B 与 nivolumab 获益相关[9]。

### MSK 条件关联经 FDR 保留，但交互作用不成立

MSK 精确重叠队列包括 138 例 RCC 和 51 个 OS 事件。73 例中位数以上 TMB 的 ccRCC 中，16 例在任一 HLA-I 位点纯合。调整年龄、性别和 ICI 方案后，任一位点纯合与较短 OS 相关（HR=4.83，95%CI 1.89-12.33，名义 P=0.00098，BH-FDR=0.0196；Figure 4A,B），比例风险检验 P=0.628，逐例删除 HR 范围为 2.94-4.05。

但交互作用证据不足。二分类 TMB-high×任一 HLA-I 纯合交互项 HR=5.44（95%CI 1.19-24.90），名义 P=0.029，BH-FDR=0.293；10 个交互项均未通过校正（Figure 4C）。最强条件结果是任一位点纯合而非发现队列的 HLA-A 特异定义，且终点为 OS 而非 PFS。因此，MSK 结果是该数据集内稳定的后验条件关联，而不是对发现状态的直接复现。

![Figure 4](figures/Figure4.png)

### 完整肽事件审计支持末端兼容性，而不支持 partner rescue

完整配对宇宙包括 63 个事件。partner allele 在 29/63 个事件中同时优于 HLA-A*02:01，占 46.0%，与 50% 无差异（双侧二项 P=0.615），故不存在全局 partner advantage。重新评分方向与归档 best-allele 标签在 63 个事件中全部一致，表明同一预测系统内的数值分配可复现。

但分配与 C 端残基强相关（Figure 5）：24/24 个 basic tail 和 3/3 个 aromatic tail 为 partner win，而 33 个 aliphatic hydrophobic tail 中仅 1 个为 partner win，其他类别为 1/3。basic 富集和 aliphatic hydrophobic 缺失在校正后仍显著（FDR=1.52×10^-12 与 7.58×10^-14）。该模式与等位基因特异 anchor preference 一致，但由于没有患者级质谱或 T 细胞实验，它仍是计算预测而非自然呈递或免疫识别证据。

![Figure 5](figures/Figure5.png)

PXD017149 中 basic C 端残基在全部配体和肿瘤特异配体中的比例分别为 17.6% 与 18.5%（Supplementary Figure S5）。这只能证明 RCC 自然 HLA 配体中可见 basic C 端，不能验证本研究的突变肽、HLA 配对或临床关联。

## 讨论

本研究不支持原先“高突变/APOBEC 底物加 HLA-A=0 构成 TKI-ICI 耐药瓶颈”的叙事。官方数据重建后，最稳妥的发现是：在 24 例队列中，3 例 HLA-A 纯合患者同时具有较高突变数和较早 PFS 事件；但突变数与 HLA-A 纯合无法在该样本中被拆分。复合标签并没有增加超出 HLA-A 纯合的独立信息，极小的 log-rank P 值不能替代可识别的交互作用。

这一边界需要正面陈述。泛癌研究已报道 HLA-I 纯合性和 HED 与 ICI 疗效的关系[11,12]，Lee 原文也已经在本发现队列和一个独立 MSK 队列报告了 mean HED 的作用[1]。本研究的增量不在于“首次发现 HLA”，而在于位点特异重建、不可识别性的直接展示、外部阴性结果的完整纳入，以及对循环筛选肽事件的纠正。

外部数据不支持可迁移的 RCC 突变-HLA 生物标志物。Miao/DFCI 的 mutation-only 不良 PFS 信号和 MSK 高 TMB 条件下的任一位点纯合不良 OS 信号值得继续验证，但 Braun、JAVELIN、Chowell 和 Nature Cancer 的结果均不支持简单、统一的方向；SNiP-RCC 又提示 HED-B 而非 HLA-A 的关联[9]。这可能由治疗方案、终点、疾病构成、测序方法和变量定义差异造成，也提示总突变数、clonal neoantigen load 与呈递范围不应互相替代。Nature Cancer 中 IO/IO 的 clonal neoantigen 与 exceptional response、IO/VEGF 的 B cell/TLS 信号也支持这种治疗背景差异[8]。

MSK 结果需要尤其谨慎。其条件关联在 20 项检验中通过 FDR，并在逐例删除中稳定；但交互项家族未通过 FDR。高 TMB 分层内存在 HLA 效应不等于 TMB 改变了 HLA 效应。该队列的治疗、OS 终点和 panel TMB 均与发现队列不同，故其仅支持进一步检验，不能称作直接验证。

完整肽审计进一步改变了机制解释。若先按 partner-only 选出 29 个事件，再声明 29/29 partner better，结论必然受到循环选择影响。纳入全部 63 个事件后，不存在全局 partner advantage，保留下来的是末端残基依赖的预测互补性。该结果可解释为什么标量化的“呈递救援”叙事并不可靠，但不能替代免疫肽组或 T 细胞功能验证。

本研究的局限性包括发现队列小、关键状态后验定义、原研究已报道 HLA 多样性、外部数据的治疗与终点不统一、MSK 结果为条件而非交互证据、肽预测依赖同一模型和输入管线且嵌套于 7 例患者，以及 PXD017149 仅为未链接的背景。研究不包含 RNA、clonality、空间免疫、患者级免疫肽组或 T 细胞读出，也无法建立治疗预测作用。后续研究应在独立、预先限定的 RCC ICI-TKI 队列中同时采集 WES、HLA、RNA、clonality、疗效和功能性抗原呈递数据。

## 结论

不同 RCC 免疫治疗队列中，突变负荷和 HLA-I 纯合性未形成一致的临床生物标志物。Lee 队列的三例观察在记录层面真实，但与 HLA-A 纯合完全重合；MSK 高 TMB 条件关联缺乏经校正的交互支持。完整肽配对审计支持末端残基依赖的等位基因互补性，而非全局 partner rescue。本文提供的是透明、可检验的证据边界与分析框架，不是可用于治疗选择的分类器。

## 数据、代码和声明

发现队列数据来自 PRJNA610643 及 Lee 补充材料；Miao/DFCI 数据来自 cBioPortal `ccrcc_dfci_2019`；MSK 数据来自 `tmb_mskcc_2018` 与 `msk_impact_50k_2026`，使用受原始许可约束；PXD017149 来自 ProteomeXchange/PRIDE。当前投稿包已附代码、运行顺序和允许发布的聚合 source data；正式投稿前仍需建立长期公开代码仓库并补入 DOI，且不得上传修改后的患者级 MSK 数据。伦理、作者贡献、经费、利益冲突及 AI 使用声明必须由作者确认，不得推测填写。

## 图注与参考文献

图注与 23 条参考文献与英文稿完全对应；投稿时应以英文稿 reference list 作为统一版本，并按目标期刊格式输出。
