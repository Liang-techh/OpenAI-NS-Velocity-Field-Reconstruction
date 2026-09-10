# 本次优化包的使用方法

本包包含已修改并测试的源代码，不是“待执行提示词”。本次没有创建远端提交、
没有推送 GitHub，也没有创建 PR。补丁针对以下已读取并按 Git 树哈希核验的版本：

```text
b01aaeebc6ec8a39f6692109b6a3c81b92cd0f32
```

保留了最初检查后新增的 5 次远端提交，包括来源清单、来源测试和分组 CI；
数值优化替代了旧的梯形积分实现，不再依赖 trapz/trapezoid 的版本差异。

## 使用已有仓库

先保留本地未提交的工作。在一个干净工作区中，将下载的补丁放到仓库外，
例如当前仓库的上一级目录。然后执行：

```bash
git switch -c improve/ns-reconstruction-numerics
git apply --check ../ns-reconstruction-improvements.patch
git apply ../ns-reconstruction-improvements.patch
python -m pip install -e ".[dev]"
python -m pytest -q -W error
python -m openai_ns_reconstruction verify
```

只有 `git apply --check` 成功才执行下一条应用命令。如果仓库在上述基准之后
又发生冲突性改动，应按冲突逐项合并，不要覆盖目录、强制重置或强制推送。
审阅代码、通过测试后，可以正常提交与推送；本包不会自动执行这些远端写入。

## 使用完整代码压缩包

解压完整源代码包，进入包含 `pyproject.toml` 的目录，执行：

```bash
python -m pip install -e ".[dev]"
python -m pytest -q -W error
python -m openai_ns_reconstruction demo --output toy-demo --grid-size 41
```

示例数据是明确标记为 **toy-not-openai** 的演示场。输出目录必须是新目录，
不会覆盖已有文件。安装时需要可用的 Python/NumPy/pytest；wheel 本身不捆绑
这些第三方依赖。

## 实际进展与验收边界

本次在 Python 3.13.5、NumPy 2.3.5 下运行：优化前初始版本 8 项测试通过；
整合并行改动后的远端基准 12 项通过；本优化版 **207 项通过、0 项失败**。
同时检查了安装、独立 wheel 导入、命令行和数值诊断。其他 Python/NumPy
组合已写入 CI 配置，但本次没有在远端执行这些 CI。

已补齐的关键部分是近奇点求根、直接 tau 接口、轴线与压力处理、径向积分、
背景场截断导数、局部化乘积项、附录 A 的逐点矩方程原语、独立残差测试和
可机器读取的验收状态。详细数据见 `reports/delivery.json`。

**尚未完成论文原速度场的全量实例化**：缺少实际主剖面、完整递归系数、
振荡波与应力实现、残差改进序列，以及穿过 t=1 的光滑外力构造。下一优先级
是 `docs/NEXT_TASKS.md` 中的 Theorem 4.6 主剖面构造，不是继续拟合 toy 场的指数。

```bash
python -m openai_ns_reconstruction verify --require-paper-exact
```

该命令当前应返回 **2**，即拒绝将部分实现判成完整论文复刻。这不是测试失败，
而是明确的完成度门槛。数值测试不能替代无限阶估计或 Lean 证明。
