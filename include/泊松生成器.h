/**
 * \file 泊松生成器.h
 * \brief
 *
 * 泊松生成器
 *
 * \version 1.6.1
 * \date 16/02/2024
 * \author Sergey Kosarevsky, 2014-2024
 * \author support@linderdaum.com   http://www.linderdaum.com   http://blog.linderdaum.com
 */

/*
   使用示例:

      #define POISSON_PROGRESS_INDICATOR 1
      #include "泊松生成器.h"
      ...
      泊松生成器::DefaultPRNG PRNG;
      const auto Points = 泊松生成器::生成泊松点集( 点数量, PRNG );
      ...
      const auto Points = 泊松生成器::生成Vogel点集( 点数量 );
*/

// 任意维度的快速泊松盘采样
// http://people.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf

// 实现基于 http://devmag.org.za/2009/05/03/poisson-disk-sampling/

/* 版本历史:
 *		1.6.1   Feb 16, 2024    使用 .clang-format 重新格式化
 *		1.6     May 29, 2023    添加 generateHammersleyPoints() 生成 Hammersley 点
 *		1.5     Mar 26, 2022    添加 generateJitteredGridPoints() 生成抖动网格点
 *		1.4.1   Dec 12, 2021		使用更快速和轻量级的 LCG 替换默认的 Mersenne Twister 和 <random>
 *		1.4     Dec  5, 2021		添加 generateVogelPoints() 生成 Vogel 盘点
 *		1.3     Mar 14, 2021		修复 bug: !isCircle 模式下的点数, 不正确的循环边界
 *		1.2     Dec 28, 2019		修复 bug; 更一致的进度指示器; 新的命令行选项在演示应用中
 *		1.1.6   Dec  7, 2019		移除重复的种子初始化; 修复警告
 *		1.1.5   Jun 16, 2019		类内初始化; 默认构造函数; 命名, 更短的代码
 *		1.1.4   Oct 19, 2016		POISSON_PROGRESS_INDICATOR 可以在头文件外部定义, 默认禁用
 *		1.1.3a  Jun  9, 2016		更新 DefaultPRNG 的构造函数
 *		1.1.3   Mar 10, 2016		头文件库, 无全局可变状态
 *		1.1.2   Apr  9, 2015		输出包含 XY 坐标的文本文件
 *		1.1.1   May 23, 2014		初始化 PRNG 种子, 修复未初始化的字段
 *		1.1     May  7, 2014		支持密度图
 *		1.0     May  6, 2014
 */

#include <stdint.h>
#include <vector>

namespace 泊松生成器 {

const char* Version = "1.6.1 (16/02/2024)";

class DefaultPRNG {
 public:
  DefaultPRNG() = default;
  explicit DefaultPRNG(unsigned int seed) : seed_(seed) {}
  inline float randomFloat() {
    seed_ *= 521167;
    uint32_t a = (seed_ & 0x007fffff) | 0x40000000;
    // remap to 0..1
    return 0.5f * (*((float*)&a) - 2.0f);
  }
  inline uint32_t randomInt(uint32_t maxInt) {
    return uint32_t(randomFloat() * maxInt);
  }
  inline uint32_t getSeed() const {
    return seed_;
  }

 private:
  uint32_t seed_ = 7133167;
};

struct 点 {
  点() = default;
  点(float X, float Y) : x(X), y(Y), 是有效的(true) {}
  float x = 0.0f;
  float y = 0.0f;
  bool 是有效的 = false;
  //
  bool 要是在矩形内() const {
    return x >= 0 && y >= 0 && x <= 1 && y <= 1;
  }
  //
  bool 要是在圆形内() const {
    const float fx = x - 0.5f;
    const float fy = y - 0.5f;
    return (fx * fx + fy * fy) <= 0.25f;
  }
  点& operator+(const 点& p) {
    x += p.x;
    y += p.y;
    return *this;
  }
  点& operator-(const 点& p) {
    x -= p.x;
    y -= p.y;
    return *this;
  }
};

struct 网格点 {
  网格点() = delete;
  网格点(int X, int Y) : x(X), y(Y) {}
  int x;
  int y;
};

float 获取距离(const 点& 起点, const 点& 终点) {
  return sqrt((起点.x - 终点.x) * (起点.x - 终点.x) + (起点.y - 终点.y) * (起点.y - 终点.y));
}

网格点 图像到网格(const 点& P, float 单格) {
  return 网格点((int)(P.x / 单格), (int)(P.y / 单格));
}

struct 网格 {
  网格(int 宽, int 高, float 单格) : 宽_(宽), 高_(高), 单格_(单格) {
    网格_.resize(高_);
    for (auto i = 网格_.begin(); i != 网格_.end(); i++) {
      i->resize(宽);
    }
  }
  void 要插入(const 点& 此点) {
    const 网格点 g = 图像到网格(此点, 单格_);
    网格_[g.x][g.y] = 此点;
  }
  bool 要是在邻近区域内(const 点& 此点, float 最小距离, float 单格尺寸) {
    const 网格点 g = 图像到网格(此点, 单格尺寸);

    // 查找邻近点的相邻单元格数量”。
    const int D = 5;

    // 扫描网格中点的邻域”。
    for (int i = g.x - D; i <= g.x + D; i++) {
      for (int j = g.y - D; j <= g.y + D; j++) {
        if (i >= 0 && i < 宽_ && j >= 0 && j < 高_) {
          const 点 P = 网格_[i][j];

          if (P.是有效的 && 获取距离(P, 此点) < 最小距离)
            return true;
        }
      }
    }

    return false;
  }

 private:
  int 宽_;
  int 高_;
  float 单格_;
  std::vector<std::vector<点>> 网格_;
};

template<typename PRNG>
点 随机取出(std::vector<点>& 点集, PRNG& 随机数生成器) {
  const int 索引 = 随机数生成器.randomInt(static_cast<int>(点集.size()) - 1);
  const 点 p = 点集[索引];
  点集.erase(点集.begin() + 索引);
  return p;
}

template<typename PRNG>
点 在周围生成随机点(const 点& 中心点, float 最小距离, PRNG& 随机数生成器) {
  // 从非均匀分布开始
  const float R1 = 随机数生成器.randomFloat();
  const float R2 = 随机数生成器.randomFloat();

  // 半径应在 最小距离 和 2 * 最小距离 之间
  const float 半径 = 最小距离 * (R1 + 1.0f);

  // 随机角度
  const float 角度 = 2 * 3.141592653589f * R2;

  // 新点围绕点 (x, y) 生成
  const float x = 中心点.x + 半径 * cos(角度);
  const float y = 中心点.y + 半径 * sin(角度);

  return 点(x, y);
}

/**
   返回生成的点集

   新增点数量 - 详细信息请参阅 bridson-siggraph07-poissondisk.pdf（值 'k'）
   是圆形  - 填充圆形则为 'true'，填充矩形则为 'false'
   最小距离 - 最小距离估计器，使用负值表示默认值
**/
template<typename PRNG = DefaultPRNG>
std::vector<点> 生成泊松点集(uint32_t 点数量, PRNG& 随机数生成器, bool 是圆形 = true, uint32_t 新增点数量 = 30, float 最小距离 = -1.0f) {
  点数量 *= 2;

  // 如果我们想要生成泊松方形形状，由于形状面积减少，将估计的点数乘以 PI/4
  if (!是圆形) {
    const double Pi_4 = 0.785398163397448309616; // PI/4
    点数量 = static_cast<int>(Pi_4 * 点数量);
  }

  if (最小距离 < 0.0f) {
    最小距离 = sqrt(float(点数量)) / float(点数量);
  }

  std::vector<点> 采样点集;
  std::vector<点> 待处理列表;

  if (!点数量)
    return 采样点集;

  // 创建网格
  const float 单格尺寸 = 最小距离 / sqrt(2.0f);

  const int 网格宽 = (int)ceil(1.0f / 单格尺寸);
  const int 网格高 = (int)ceil(1.0f / 单格尺寸);

  网格 网格值(网格宽, 网格高, 单格尺寸);

  点 首个点;
  do {
    首个点 = 点(随机数生成器.randomFloat(), 随机数生成器.randomFloat());
  } while (!(是圆形 ? 首个点.要是在圆形内() : 首个点.要是在矩形内()));

  // 更新容器
  待处理列表.push_back(首个点);
  采样点集.push_back(首个点);
  网格值.要插入(首个点);

#if POISSON_PROGRESS_INDICATOR
  size_t progress = 0;
#endif

  // 为队列中的每个点生成新点。
  while (!待处理列表.empty() && 采样点集.size() <= 点数量) {
#if POISSON_PROGRESS_INDICATOR
    // a progress indicator, kind of
    if ((samplePoints.size()) % 1000 == 0) {
      const size_t newProgress = 200 * (samplePoints.size() + processList.size()) / numPoints;
      if (newProgress != progress) {
        progress = newProgress;
        std::cout << ".";
      }
    }
#endif // POISSON_PROGRESS_INDICATOR

    const 点 当前点 = 随机取出<PRNG>(待处理列表, 随机数生成器);

    for (uint32_t i = 0; i < 新增点数量; i++) {
      const 点 新点 = 在周围生成随机点(当前点, 最小距离, 随机数生成器);
      const bool 是可放置点 = 是圆形 ? 新点.要是在圆形内() : 新点.要是在矩形内();

      if (是可放置点 && !网格值.要是在邻近区域内(新点, 最小距离, 单格尺寸)) {
        待处理列表.push_back(新点);
        采样点集.push_back(新点);
        网格值.要插入(新点);
        continue;
      }
    }
  }

#if POISSON_PROGRESS_INDICATOR
  std::cout << std::endl << std::endl;
#endif // POISSON_PROGRESS_INDICATOR

  return 采样点集;
}

点 采样Vogel盘(uint32_t 索引, uint32_t 点数量, float 角度) {
  const float 黄金角 = 2.4f;

  const float 半径 = sqrtf(float(索引) + 0.5f) / sqrtf(float(点数量));
  const float 新角度 = 索引 * 黄金角 + 角度;

  return 点(半径 * cosf(新角度), 半径 * sinf(新角度));
}

/**
  返回生成的点集
**/
std::vector<点> 生成Vogel点集(uint32_t 点数量, bool 是圆形 = true, float 角度 = 0.0f, 点 中心点 = 点(0.5f, 0.5f)) {
  std::vector<点> 采样点集;

  采样点集.reserve(点数量);

  const uint32_t 采样数 = 是圆形 ? 4 * 点数量 : 点数量;

  for (uint32_t i = 0; i != 点数量; i++) {
    const 点 新点 = 采样Vogel盘(i, 采样数, 角度 * 3.141592653f / 180.0f) + 中心点;
    采样点集.push_back(新点);
  }

  return 采样点集;
}

/**
  返回生成的点向量

  泊松盘 VS 抖动网格 https://www.redblobgames.com/x/1830-jittered-grid/
**/
template<typename PRNG = DefaultPRNG>
std::vector<点> 生成抖动网格点集(uint32_t 点数量,
                                 PRNG& 随机数生成器,
                                 bool 是圆形 = false,
                                 float 抖动半径 = 0.004f,
                                 点 中心点 = 点(0.5f, 0.5f)) {
  std::vector<点> 采样点集;

  采样点集.reserve(点数量);

  const uint32_t 网格尺寸 = uint32_t(sqrt(点数量));

  for (uint32_t x = 0; x != 网格尺寸; x++) {
    for (uint32_t y = 0; y != 网格尺寸; y++) {
      点 新点;
      do {
        const 点 偏移点 = 在周围生成随机点(点(0, 0), 抖动半径, 随机数生成器) - 中心点 + 点(0.5f, 0.5f);
        新点 = 点(float(x) / 网格尺寸, float(y) / 网格尺寸) + 偏移点;
        // 生成一个新点，直到它在边界内
      } while (!新点.要是在矩形内());

      if (是圆形)
        if (!新点.要是在圆形内())
          continue;

      采样点集.push_back(新点);
    }
  }

  return 采样点集;
}

namespace {

// http://holger.dammertz.org/stuff/notes_HammersleyOnHemisphere.html
float radicalInverse_VdC(uint32_t bits) {
  bits = (bits << 16u) | (bits >> 16u);
  bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
  bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
  bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
  bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
  return float(float(bits) * 2.3283064365386963e-10); // / 0x100000000
}

点 hammersley2d(uint32_t i, uint32_t N) {
  return 点(float(i) / float(N), radicalInverse_VdC(i));
}

} // namespace

/**
  返回生成的点集
**/
std::vector<点> 生成Hammersley点集(uint32_t 点数量) {
  std::vector<点> 采样点集;

  采样点集.reserve(点数量);

  const uint32_t gridSize = uint32_t(sqrt(点数量));

  for (uint32_t i = 0; i != 点数量; i++) {
    点 p = hammersley2d(i, 点数量);

    采样点集.push_back(p);
  }
  return 采样点集;
}

} // namespace 泊松生成器
