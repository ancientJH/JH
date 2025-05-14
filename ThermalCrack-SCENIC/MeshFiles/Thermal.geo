// Gmsh 脚本：方块，底部通过尺寸场加密网格 (单一表面)

// 特征长度定义
ell = 7.5e-5;
lc1 = ell/3;   // 细网格尺寸 (目标加密区域)
lc2 = ell;     // 中等网格尺寸 (此脚本中未使用)
lc3 = ell*6;   // 粗网格尺寸 (背景或远离加密区域)

// 几何尺寸
L = 5e-3;
H = 5e-3;
H_quarter = H / 4; // 加密区域的界限高度

// 定义点，并在定义时直接附加期望的网格特征长度
// Point(tag) = {x, y, z, characteristic_length};
Point(1) = {0, 0, 0, lc1};           // 左下角点 P1
Point(2) = {L, 0, 0, lc1};           // 右下角点 P2
Point(3) = {L, H, 0, lc3};           // 右上角点 P3
Point(4) = {0, H, 0, lc3};           // 左上角点 P4

// 定义构成外部边界的线
Line(1) = {1, 2}; // 底边
Line(2) = {2, 3}; // 右边
Line(3) = {3, 4}; // 顶边
Line(4) = {4, 1}; // 左边

// 定义单一的平面覆盖整个区域
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// 定义物理组
Physical Curve("top", 10) = {3};
Physical Curve("bottom", 11) = {1};
Physical Curve("left", 12) = {4};
Physical Curve("right", 13) = {2};
Physical Surface("domain", 20) = {1}; // 整个计算域

// --- 使用背景网格尺寸场来控制加密 ---

// 1. 定义一个 Box 场，用于指定底部 H/4 区域的网格尺寸
Field[1] = Box;
Field[1].VIn = lc1;       // Box 内部的网格尺寸 (细)
Field[1].VOut = lc3;      // Box 外部的网格尺寸 (粗) - 这个值会与点定义的尺寸通过Min场结合
Field[1].XMin = 0;
Field[1].XMax = L;
Field[1].YMin = 0;
Field[1].YMax = H_quarter; // Box 的上边界即为虚拟分界线
Field[1].ZMin = -L/2;     // Z范围需要覆盖2D平面
Field[1].ZMax = L/2;
// Field[1].Thickness = H; // 可选，定义过渡区域的厚度，如果需要更平滑的过渡

// 2. （可选但推荐）定义一个 Threshold 场，可以基于点定义的特征长度来影响网格
//    这里我们主要依赖 Box 场的 VOut 和点自身的特征长度，
//    并通过 Min 场来确保点定义的 lc1 和 lc3 被尊重。
//    如果只用 Box 场，VOut 会决定 Box 外的尺寸。
//    为了让点定义的 lc1 和 lc3 也能在 Box 区域外及边界处起作用，
//    我们可以直接将 Box 场作为背景场，Gmsh 会自动插值。
//    或者，更精确地，使用 Min 场结合一个从点插值的场。

// 为了简单起见，我们先尝试仅用 Box 场，并依赖点上定义的特征长度进行插值。
// Gmsh 在生成网格时，会考虑点上定义的特征长度，并从这些点向内部插值。
// Box 场会提供一个额外的约束。

// 将 Box 场设置为背景场
Background Field = 1;

// 或者，使用 Min 场来组合多个影响因素（如果需要更复杂的控制）
// Field[2] = Min;
// Field[2].FieldsList = {1}; // 包含 Box 场
// Background Field = 2;


Coherence; // 确保几何定义的一致性

// 网格生成命令
Mesh 2; // 生成二维网格
