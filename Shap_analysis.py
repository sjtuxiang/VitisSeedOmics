import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import seaborn as sns
import shap
import networkx as nx
from mpl_toolkits.axes_grid1 import make_axes_locatable
import statsmodels.api as sm
import warnings
warnings.filterwarnings("ignore")

# --- 模型与预处理库 ---
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.feature_selection import f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier

# ==========================================
# 1. 全局配置与参数设置
# ==========================================
CONFIG = {
    "data_file": "data.txt",         
    "output_dir": "ML_Results",    
    "dpi": 300,
    "formats": ["pdf", "png"],    
    
    "font_family_en": "sans-serif",  
    "label_fontsize": 14,
    "title_fontsize": 16,
    "tick_fontsize": 14,
    
    # === 决定进入最终分析与所有图表展示的核心物质数量 ===
    "top_k_features_to_retain": 50,  
    "fig5_max_features": 20,         
    "fig7_max_features": 20,         
    "fig8_max_display": 20,          
    
    "fig1_kde_linewidth": 1.5,
    "fig1_class0_color": "#6b6aa3",  
    "fig1_class1_color": "#54748a",  
    "fig2_scatter_s": 12,
    "fig2_bar_color": "#6b6aa3",
    "fig3_scatter_s": 20,
    "fig3_curve_color": "#54748a",
    "fig3_curve_linewidth": 1.5,
    "fig4_main_color": "#54748a",
    "fig4_inter_color": "#ad4447",
    "fig7_node_color": "#6b6aa3",    
    "fig7_edge_color": "#ad4447",    
    "fig8_figsize": (12, 10),
    "fig10_line_color_q1": "#54748a", 
    "fig10_line_color_med": "#6b6aa3",
    "fig10_line_color_q3": "#ad4447", 
    
    "shap_cmap": mcolors.LinearSegmentedColormap.from_list("custom_cmap", [ "#54748a", "#FFFFFF", "#ad4447"]),
    "random_state": 42
}

import logging
# 关闭 Matplotlib 的字体查找警告
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['pdf.fonttype'] = 42

# ==========================================
# 2. 辅助函数
# ==========================================
def save_figure(fig, filename_base, subfolder=None):
    save_dir = CONFIG["output_dir"]
    if subfolder:
        save_dir = os.path.join(save_dir, subfolder)
    os.makedirs(save_dir, exist_ok=True)
    safe_filename_base = str(filename_base).replace('/', '_').replace('\\', '_')
    for fmt in CONFIG["formats"]:
        filepath = os.path.join(save_dir, f"{safe_filename_base}.{fmt}")
        fig.savefig(filepath, dpi=CONFIG["dpi"], format=fmt, bbox_inches='tight')
    plt.close(fig)

# ==========================================
# 3. 核心分析类构建
# ==========================================
class XGBoostXAIAnalyzer:
    def __init__(self):
        self.random_state = CONFIG["random_state"]
        self.model = None
        self.pipeline = None
        self.explainer = None
        self.shap_values_obj = None
        self.shap_interaction_values = None
        self.X_train_processed = None

    def prepare_data(self):
        print("正在读取代谢组学数据表...")
        df = pd.read_csv(CONFIG["data_file"], sep='\t')
        
        # 精准提取仅包含 Vv_ 和 Vh_ 的样本列，绕开所有文字注释列
        sample_cols = [col for col in df.columns if 'Vv_' in str(col) or 'Vh_' in str(col)]
        
        # 提取第一列(Compounds)作为特征名
        raw_names = df.iloc[:, 0].values.tolist()
        self.feature_names = [str(name).replace('[', '(').replace(']', ')').replace('<', '_') for name in raw_names]
        
        # 转置数据，并强制转换为 float 浮点型
        self.X = pd.DataFrame(df[sample_cols].values.T, columns=self.feature_names).astype(float)
        
        # 动态生成标签：只要列名中包含 'Vh' 就视为 1，包含 'Vv' 视为 0
        labels = [1 if 'Vh' in str(col) else 0 for col in sample_cols]
        self.y = pd.Series(labels)
        
        # 剔除零方差特征（防止缩放和方差分析时报错）
        variances = self.X.var()
        non_zero_var_cols = variances[variances > 0].index
        self.X = self.X[non_zero_var_cols]
        self.feature_names = self.X.columns.tolist()
        
        self.X_train, self.y_train = self.X, self.y
        print(f"数据读取完毕。总样本数: {self.X.shape[0]}, 有效特征数量(剔除零方差): {self.X.shape[1]}")

    def ensemble_feature_screening(self):
        print("🚀 正在使用 ANOVA单变量 + ML模型 进行加权融合特征筛选...")
        
        # -----------------------------------------------------------------
        # 新增核心功能：计算原始数据的均值和 Log2 Fold Change (Log2FC)
        # -----------------------------------------------------------------
        # label=1 是 Vh, label=0 是 Vv
        mean_Vh = self.X_train[self.y_train == 1].mean(axis=0)
        mean_Vv = self.X_train[self.y_train == 0].mean(axis=0)
        
        # 计算 Log2FC，加一个小常数(epsilon)防止除以0或取对数报错
        epsilon = 1e-9
        log2fc = np.log2((mean_Vh + epsilon) / (mean_Vv + epsilon))
        # -----------------------------------------------------------------

        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(self.X_train)

        importances = pd.DataFrame(index=self.feature_names)
        
        # 将均值和 Log2FC 记录到表格中
        importances["Mean_Vv"] = mean_Vv.values
        importances["Mean_Vh"] = mean_Vh.values
        importances["Log2FC(Vh/Vv)"] = log2fc.values

        minmax = MinMaxScaler()

        # 引入单变量统计检验 (ANOVA F-score)
        f_vals, p_vals = f_classif(X_scaled, self.y_train)
        f_vals = np.nan_to_num(f_vals) # 防止除零错误产生 NaN
        importances["Univariate_F"] = minmax.fit_transform(f_vals.reshape(-1, 1)).flatten()
        importances["P_Value"] = p_vals # 顺便记录 p-value 供参考

        models = {
            "LR": LogisticRegression(penalty='l2', class_weight='balanced', random_state=self.random_state),
            "SVM": SVC(kernel='linear', class_weight='balanced', random_state=self.random_state),
            "RF": RandomForestClassifier(n_estimators=300, max_depth=3, max_features='sqrt', class_weight='balanced', random_state=self.random_state),
            "XGB": XGBClassifier(n_estimators=100, max_depth=2, learning_rate=0.03, subsample=0.7, colsample_bytree=0.5, eval_metric='logloss', use_label_encoder=False, random_state=self.random_state),
            "LGBM": LGBMClassifier(n_estimators=100, max_depth=3, colsample_bytree=0.3, class_weight='balanced', min_child_samples=2, min_split_gain=0.0, verbose=-1, random_state=self.random_state)
        }

        for name, model in models.items():
            model.fit(X_scaled, self.y_train)
            if name in ["LR", "SVM"]:
                imp = np.abs(model.coef_[0])
            else:
                imp = model.feature_importances_
            
            importances[name] = minmax.fit_transform(imp.reshape(-1, 1)).flatten()

        # 加权共识得分
        importances['Consensus_Score'] = (
            importances["Univariate_F"] * 0.75 + 
            importances["LR"] * 0.05 + 
            importances["SVM"] * 0.05 + 
            importances["RF"] * 0.05 + 
            importances["XGB"] * 0.05 + 
            importances["LGBM"] * 0.05
        )
        
        # -----------------------------------------------------------------
        # 新增核心功能：根据方向进行过滤，只保留 Vh 显著高于 Vv 的物质
        # -----------------------------------------------------------------
        # Log2FC > 0 代表 Vh > Vv。
        up_in_Vh_mask = importances["Log2FC(Vh/Vv)"] > 1
        
        # 分离出两张表：全量表，以及“仅Vh高”的表
        importances_all = importances.sort_values(by='Consensus_Score', ascending=False)
        importances_up_in_Vh = importances[up_in_Vh_mask].sort_values(by='Consensus_Score', ascending=False)

        os.makedirs(CONFIG["output_dir"], exist_ok=True)
        # 保存全量特征供参考
        importances_all.to_csv(os.path.join(CONFIG["output_dir"], "MultiModel_Feature_Ranking_ALL.csv"))
        # 保存过滤后的特征
        importances_up_in_Vh.to_csv(os.path.join(CONFIG["output_dir"], "MultiModel_Feature_Ranking_Vh_Upregulated.csv"))
        
        print(f"✅ 融合特征筛选完成！过滤前物质总数: {len(importances_all)}，其中 Vh > Vv 的物质数: {len(importances_up_in_Vh)}")
        
        # 强制只从 Vh > Vv 的物质中提取前 Top K 传入后续的 SHAP 分析
        top_features = importances_up_in_Vh.index[:CONFIG["top_k_features_to_retain"]].tolist()
        
        if len(top_features) == 0:
            raise ValueError("没有找到 Vh 含量大于 Vv 的特征，请检查数据！")

        self.feature_names = top_features
        self.X_train = self.X_train[top_features]
        print(f"✅ 已精准截取 Top {len(top_features)} 个 [Vh高表达] 的高价值物质用于后续 SHAP 分析。")

    def train_final_model(self):
        print("正在训练最终解释模型 (XGBoost) - 已启用抗特征垄断机制...")
        xgb_clf = XGBClassifier(
            n_estimators=300,       
            max_depth=2,             
            learning_rate=0.01,      
            subsample=0.7,           
            colsample_bytree=0.3,    
            colsample_bylevel=0.5,   
            reg_alpha=1.5,           
            reg_lambda=3.0,          
            eval_metric='logloss',
            use_label_encoder=False,
            random_state=self.random_state
        )
        
        self.pipeline = Pipeline([
            ('scaler', StandardScaler()),
            ('classifier', xgb_clf) 
        ])
        
        self.pipeline.fit(self.X_train, self.y_train)
        self.model = self.pipeline.named_steps['classifier']

    def calculate_shap(self):
        print("正在计算 SHAP 及交互值...")
        scaler = self.pipeline.named_steps['scaler']
        X_scaled = scaler.transform(self.X_train)
        self.X_train_processed = pd.DataFrame(X_scaled, columns=self.feature_names)
        
        self.explainer = shap.TreeExplainer(self.model)
        self.shap_values_obj = self.explainer(self.X_train_processed)
        self.shap_interaction_values = self.explainer.shap_interaction_values(self.X_train_processed)
        print("SHAP 计算完成")

    # ================= 图表绘制部分 =================
    def plot_figure_1(self):
        print("1/7 正在绘制图1：模型分类性能图...")
        preds_proba = self.pipeline.predict_proba(self.X_train)[:, 1]
        y_true = self.y_train.values  
        
        fig, ax = plt.subplots(figsize=(7, 5))
        sns.kdeplot(preds_proba[y_true==0], fill=True, color=CONFIG["fig1_class0_color"], ax=ax)
        sns.kdeplot(preds_proba[y_true==1], fill=True, color=CONFIG["fig1_class1_color"], ax=ax)
        sns.scatterplot(x=preds_proba, y=[0.02]*len(preds_proba), hue=y_true, 
                        palette={0:CONFIG["fig1_class0_color"], 1:CONFIG["fig1_class1_color"]}, 
                        s=120, zorder=5, legend=False, ax=ax, edgecolor='white')
        
        patch1 = mlines.Line2D([], [], color=CONFIG["fig1_class0_color"], marker='o', linestyle='None', markersize=10, label='Vv_EL (Class 0)')
        patch2 = mlines.Line2D([], [], color=CONFIG["fig1_class1_color"], marker='o', linestyle='None', markersize=10, label='Vh_EL (Class 1)')
        ax.legend(handles=[patch1, patch2], frameon=False, fontsize=CONFIG["label_fontsize"])
        ax.set_title("Fig 1: Classification Probability\n(Vv vs Vh)", fontsize=CONFIG["title_fontsize"], fontweight='bold')
        ax.set_xlabel("Predicted Probability of being 'Vh_EL'", fontsize=CONFIG["label_fontsize"])
        ax.set_ylabel("Density", fontsize=CONFIG["label_fontsize"])
        save_figure(fig, "Fig1_Classification_Performance")

    def plot_figure_2(self):
        print("2/7 正在绘制图2：全局特征贡献度分析图...")
        fig, ax1 = plt.subplots(figsize=(10, 12))
        mean_abs_shap = np.abs(self.shap_values_obj.values).mean(axis=0)
        sort_inds = np.argsort(mean_abs_shap)[-20:]
        sorted_features = [self.feature_names[i] for i in sort_inds]
        sorted_mean_shap = mean_abs_shap[sort_inds]
        total_shap_sum = np.sum(mean_abs_shap)
        shap_vals = self.shap_values_obj.values[:, sort_inds]
        feat_vals = self.X_train_processed.values[:, sort_inds]
        y_pos = np.arange(len(sorted_features))
        
        ax2 = ax1.twiny()
        ax2.barh(y_pos, sorted_mean_shap, color=CONFIG["fig2_bar_color"], align='center', alpha=0.8, height=0.6, zorder=2)
        cmap = CONFIG["shap_cmap"]
        for i in range(len(sorted_features)):
            row_shap, row_feat = shap_vals[:, i], feat_vals[:, i]
            fmin, fmax = np.min(row_feat), np.max(row_feat)
            row_feat_norm = (row_feat - fmin) / (fmax - fmin) if fmax > fmin else np.zeros_like(row_feat)
            jitter = np.random.normal(0, 0.1, size=len(row_shap))
            ax1.scatter(row_shap, np.repeat(i, len(row_shap)) + jitter, c=row_feat_norm, cmap=cmap, s=CONFIG["fig2_scatter_s"], alpha=0.8, edgecolors='none', zorder=4)
            
        max_mean_val = np.max(sorted_mean_shap)
        ax2.set_xlim(0, max_mean_val * 1.2)
        for i, v in enumerate(sorted_mean_shap):
            ax1.text(v + max_mean_val * 0.01, i, f"{(v / total_shap_sum) * 100:.2f}%", va='center', ha='left', fontsize=CONFIG["label_fontsize"]-2, transform=ax2.transData, zorder=10)
        
        ax1.set_zorder(ax2.get_zorder() + 1)
        ax1.patch.set_visible(False)
        ax1.set_yticks(y_pos)
        ax1.set_yticklabels(sorted_features, fontsize=CONFIG["tick_fontsize"])
        ax1.set_xlabel("SHAP value (impact on model output)", fontsize=CONFIG["label_fontsize"])
        ax2.set_xlabel("Mean Absolute SHAP Value", fontsize=CONFIG["label_fontsize"])
        
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="3%", pad=0.1)
        cbar = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1)), cax=cax)
        cbar.set_ticks([0, 1])
        cbar.set_ticklabels(['Low', 'High'])
        cbar.set_label('Feature value', rotation=270, labelpad=15, fontsize=CONFIG["label_fontsize"])
        save_figure(fig, "Fig2_Global_Contribution")

    def plot_figure_3(self):
        print("3/7 正在绘制图3：单特征偏依赖图...")
        mean_abs_shap = np.abs(self.shap_values_obj.values).mean(axis=0)
        sort_inds_desc = np.argsort(mean_abs_shap)[::-1][:20]
        folder_name = "Fig3_Dependence_Plots"
        for rank, feat_idx in enumerate(sort_inds_desc):
            feature_name = self.feature_names[feat_idx]
            fig, ax = plt.subplots(figsize=(6, 5))
            x_vals = self.X_train_processed[feature_name].values
            y_vals = self.shap_values_obj.values[:, feat_idx]
            scatter = ax.scatter(x_vals, y_vals, c=y_vals, cmap=CONFIG["shap_cmap"], s=CONFIG["fig3_scatter_s"], alpha=0.9, edgecolors='none', zorder=4)
            sorted_indices = np.argsort(x_vals)
            x_sorted, y_sorted = x_vals[sorted_indices], y_vals[sorted_indices]
            try:
                lowess = sm.nonparametric.lowess(y_sorted, x_sorted, frac=0.4)
                ax.plot(lowess[:, 0], lowess[:, 1], color=CONFIG["fig3_curve_color"], linewidth=CONFIG["fig3_curve_linewidth"], label="Lowess curve", zorder=5)
            except: pass
            
            ax.set_xlabel("Feature value (Z-score)", fontsize=CONFIG["label_fontsize"])
            ax.set_ylabel("SHAP Value", fontsize=CONFIG["label_fontsize"])
            ax.set_title(feature_name[:30], fontsize=CONFIG["title_fontsize"])
            ax.axhline(0, color='gray', linestyle='--', alpha=0.5, zorder=2)
            cbar = fig.colorbar(scatter, ax=ax)
            cbar.set_ticks([y_vals.min(), y_vals.max()])
            cbar.set_ticklabels(['Low', 'High'])
            cbar.set_label('SHAP value', rotation=270, labelpad=15)
            clean_feat = re.sub(r'[^\w]', '', feature_name)
            save_figure(fig, f"{rank+1:02d}_{clean_feat}_Dependence", subfolder=folder_name)

    def plot_figure_4(self):
        print("4/7 正在绘制图4：主效应与交互效应对比图...")
        num_features = min(20, len(self.feature_names))
        mean_abs_shap = np.abs(self.shap_values_obj.values).mean(axis=0)
        sort_inds = np.argsort(mean_abs_shap)[::-1][:num_features]
        main_effects = np.zeros(num_features)
        inter_effects = np.zeros(num_features)
        for i, idx in enumerate(sort_inds):
            main_effects[i] = np.abs(self.shap_interaction_values[:, idx, idx]).mean()
            mask = np.ones(len(self.feature_names), dtype=bool)
            mask[idx] = False
            inter_effects[i] = np.abs(self.shap_interaction_values[:, idx, mask]).sum(axis=1).mean()
            
        top_names = [self.feature_names[i][:20] for i in sort_inds]
        fig, ax = plt.subplots(figsize=(12, 6))
        x = np.arange(len(top_names))
        width = 0.4
        ax.bar(x - width/2, main_effects, width, label='Main effect', color=CONFIG["fig4_main_color"])
        ax.bar(x + width/2, inter_effects, width, label='Interaction', color=CONFIG["fig4_inter_color"])
        ax.set_xticks(x)
        ax.set_xticklabels(top_names, rotation=45, ha='right', fontsize=9)
        ax.legend()
        ax.set_title('Main vs Interaction Effect Magnitude', fontweight='bold')
        save_figure(fig, "Fig4_Main_vs_Interaction")

    def plot_figure_5_7(self):
        print("5&6/7 正在绘制图5/7：交互矩阵与网络图...")
        num_top = min(CONFIG["fig5_max_features"], len(self.feature_names))
        mean_abs_shap = np.abs(self.shap_values_obj.values).mean(axis=0)
        top_inds = np.argsort(mean_abs_shap)[::-1][:num_top]
        
        top_names_full = [self.feature_names[i] for i in top_inds] 
        top_names_short = [name[:20] for name in top_names_full]
        
        inter_matrix = np.zeros((num_top, num_top))
        for i in range(num_top):
            for j in range(num_top):
                if i != j:
                    inter_matrix[i, j] = np.abs(self.shap_interaction_values[:, top_inds[i], top_inds[j]]).mean()
        
        fig5 = plt.figure(figsize=(12, 10))
        sns.heatmap(pd.DataFrame(inter_matrix, index=top_names_short, columns=top_names_short), 
                    cmap=CONFIG["shap_cmap"], annot=False, cbar_kws={'label': 'Mean |SHAP Interaction|'}) # annot设为False防止过挤
        plt.title("Fig 5: Interaction Matrix", fontweight='bold')
        save_figure(fig5, "Fig5_Interaction_Matrix")
        
        fig7 = plt.figure(figsize=(10, 10))
        G = nx.Graph()
        edge_data = []
        for i in range(num_top):
            for j in range(i+1, num_top):
                edge_data.append((i, j, inter_matrix[i, j]))
                
        weights_arr = [x[2] for x in edge_data]
        threshold = np.percentile(weights_arr, 50) 
        
        for i, name_full in enumerate(top_names_full):
            main_e = np.abs(self.shap_interaction_values[:, top_inds[i], top_inds[i]]).mean()
            G.add_node(name_full, raw_size=main_e, short_name=top_names_short[i])
            
        node_sizes_raw = [G.nodes[n]['raw_size'] for n in G.nodes()]
        max_ns = max(node_sizes_raw) if max(node_sizes_raw) > 0 else 1.0
        norm_node_sizes = [600 + 1500 * (s / max_ns) for s in node_sizes_raw] 
        
        for (i, j, w) in edge_data:
            if w >= threshold and w > 1e-8: 
                G.add_edge(top_names_full[i], top_names_full[j], raw_weight=w)
                
        if len(G.edges) > 0:
            edge_weights_raw = [G[u][v]['raw_weight'] for u, v in G.edges]
            max_ew = max(edge_weights_raw) if max(edge_weights_raw) > 0 else 1.0
            min_ew = min(edge_weights_raw)
            norm_edge_widths = [1.0 + 4.0 * ((ew - min_ew) / (max_ew - min_ew + 1e-8)) for ew in edge_weights_raw]
        else:
            norm_edge_widths = []

        pos = nx.circular_layout(G)
        nx.draw_networkx_nodes(G, pos, node_size=norm_node_sizes, node_color=CONFIG["fig7_node_color"], edgecolors='white', linewidths=2, alpha=0.9)
        if len(G.edges) > 0:
            nx.draw_networkx_edges(G, pos, width=norm_edge_widths, edge_color=CONFIG["fig7_edge_color"], alpha=0.6)
            
        labels = {n: G.nodes[n]['short_name'] for n in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, font_weight='bold')
        
        plt.title("Fig 7: Interaction Network", fontweight='bold', fontsize=16)
        plt.axis('off')
        plt.margins(0.15)
        save_figure(fig7, "Fig7_Interaction_Network")

    def plot_figure_8(self):
        print("7/7 正在绘制图8：SHAP 样本热力图...")
        fig = plt.figure(figsize=CONFIG["fig8_figsize"])
        mean_abs_shap = np.abs(self.shap_values_obj.values).mean(axis=0)
        top_idx = np.argsort(mean_abs_shap)[::-1][:CONFIG["fig8_max_display"]]
        
        preds_proba = self.pipeline.predict_proba(self.X_train)[:, 1]
        sample_order = np.argsort(preds_proba)
        ordered_shap_values = self.shap_values_obj.values[sample_order]
        
        sns.heatmap(ordered_shap_values[:, top_idx].T, cmap=CONFIG["shap_cmap"], 
                    yticklabels=[self.feature_names[i][:20] for i in top_idx])
        plt.title("Fig 8: SHAP Heatmap (Samples sorted by Pred Proba)", fontweight='bold')
        plt.xlabel("Samples (Sorted by Model Probability)")
        save_figure(fig, "Fig8_SHAP_Heatmap")

    def plot_figure_10(self):
        print("正在绘制图10：二维PDP特征依赖轮廓图...")
        mean_abs_shap = np.abs(self.shap_values_obj.values).mean(axis=0)
        sort_inds = np.argsort(mean_abs_shap)[::-1]
        X_bg_median = self.X_train_processed.median().values
        grid_resolution = 50
        cmap = CONFIG["shap_cmap"]
        
        # --- 取前 20 个最重要的特征绘制核心的两两交互图 ---
        top_n_for_pdp = 20 
        top_indices = sort_inds[:top_n_for_pdp]
        
        for i in range(len(top_indices)):
            for j in range(i + 1, len(top_indices)): # i < j 避免生成重复翻转的图
                idx1 = top_indices[i]
                idx2 = top_indices[j]
                
                feat1 = self.feature_names[idx1]
                feat2 = self.feature_names[idx2]
                clean_feat1 = re.sub(r'[^\w]', '', feat1)
                clean_feat2 = re.sub(r'[^\w]', '', feat2)
                folder_name = f"Fig10_2D_PDP"
                os.makedirs(os.path.join(CONFIG["output_dir"], folder_name), exist_ok=True)
                
                x1_min, x1_max = self.X_train_processed[feat1].min(), self.X_train_processed[feat1].max()
                x1_grid = np.linspace(x1_min, x1_max, grid_resolution)
                
                x2_min, x2_max = self.X_train_processed[feat2].min(), self.X_train_processed[feat2].max()
                x2_grid = np.linspace(x2_min, x2_max, grid_resolution)
                
                X1, X2 = np.meshgrid(x1_grid, x2_grid)
                grid_points = np.tile(X_bg_median, (grid_resolution * grid_resolution, 1))
                grid_points[:, idx1] = X1.ravel()
                grid_points[:, idx2] = X2.ravel()
                
                preds = self.model.predict_proba(grid_points)[:, 1].reshape(grid_resolution, grid_resolution)
                fig, ax = plt.subplots(figsize=(7, 6))
                cf = ax.contourf(X1, X2, preds, levels=30, cmap=cmap, alpha=0.85)
                
                q1, med, q3 = np.percentile(preds, [25, 50, 75])
                p_max_val, p_min_val = preds.max(), preds.min()
                max_idx = np.unravel_index(np.argmax(preds), preds.shape)
                min_idx = np.unravel_index(np.argmin(preds), preds.shape)
                
                try:
                    ax.contour(X1, X2, preds, levels=[q1], colors=[CONFIG["fig10_line_color_q1"]], linestyles=['--'], linewidths=1.5)
                    ax.contour(X1, X2, preds, levels=[med], colors=[CONFIG["fig10_line_color_med"]], linestyles=['-'], linewidths=1.5)
                    ax.contour(X1, X2, preds, levels=[q3], colors=[CONFIG["fig10_line_color_q3"]], linestyles=['--'], linewidths=1.5)
                except ValueError: pass 
                
                ax.scatter(X1[max_idx], X2[max_idx], marker='*', color='orange', s=180, edgecolors='black', linewidth=1, zorder=5)
                ax.scatter(X1[min_idx], X2[min_idx], marker='o', color='#00BCD4', s=80, edgecolors='white', linewidth=1.5, zorder=5)
                cbar = fig.colorbar(cf, ax=ax, pad=0.03)
                cbar.set_label("Probability of 'Vh_EL'")
                cbar.ax.tick_params(labelsize=CONFIG["tick_fontsize"])
                
                l1 = mlines.Line2D([], [], color=CONFIG["fig10_line_color_q1"], linestyle='--', label=f'Q1: {q1:.3f}')
                l2 = mlines.Line2D([], [], color=CONFIG["fig10_line_color_med"], linestyle='-', label=f'Median: {med:.3f}')
                l3 = mlines.Line2D([], [], color=CONFIG["fig10_line_color_q3"], linestyle='--', label=f'Q3: {q3:.3f}')
                p_max_leg = ax.scatter([], [], marker='*', color='orange', s=100, edgecolors='black', label=f'Max: {p_max_val:.3f}')
                p_min_leg = ax.scatter([], [], marker='o', color='#00BCD4', s=60, edgecolors='white', label=f'Min: {p_min_val:.3f}')
                
                leg = ax.legend(handles=[l1, l2, l3, p_max_leg, p_min_leg], loc='best', fontsize=CONFIG["label_fontsize"] - 3, frameon=True, facecolor=(1, 1, 1, 0.85), edgecolor='gray')
                for text in leg.get_texts(): text.set_color('black')
                
                ax.set_xlabel(feat1[:30], fontsize=CONFIG["label_fontsize"])
                ax.set_ylabel(feat2[:30], fontsize=CONFIG["label_fontsize"])
                ax.set_title("Fig 10: 2D Partial Dependence", fontweight='bold', fontsize=CONFIG["title_fontsize"])
                save_figure(fig, f"{clean_feat1}_vs_{clean_feat2}", subfolder=folder_name)

    def run_all(self):
        self.prepare_data()
        self.ensemble_feature_screening()
        self.train_final_model()
        self.calculate_shap()
        
        self.plot_figure_1()
        self.plot_figure_2()
        self.plot_figure_3()
        self.plot_figure_4()
        self.plot_figure_5_7()
        self.plot_figure_8()
        self.plot_figure_10() 
        
        print(f"🎉 全部模型训练及图表绘制完成！文件已存至 {CONFIG['output_dir']} 目录！")

if __name__ == "__main__":
    analyzer = XGBoostXAIAnalyzer()
    analyzer.run_all()
