

#include "raylib.h"
#include "泊松生成器.h"

//------------------------------------------------------------------------------------
// 入口
//------------------------------------------------------------------------------------
int main(void) {
  // 初始化
  //--------------------------------------------------------------------------------------
  const int 屏幕宽 = 800;
  const int 屏幕高 = 450;

  InitWindow(屏幕宽, 屏幕高, "泊松生成器");

  SetTargetFPS(60); // 设置帧速率
  //--------------------------------------------------------------------------------------
  const int 点数量 = 100;
  泊松生成器::DefaultPRNG PRNG; // 随机数生成器
  const auto 泊松点集 = 泊松生成器::生成泊松点集(点数量, PRNG);
  const auto 抖动网格点集 = 泊松生成器::生成抖动网格点集(100, PRNG, true, 0.015f);
  const auto Vogel点集 = 泊松生成器::生成Vogel点集(点数量);
  const auto Hammersley点集 = 泊松生成器::生成Hammersley点集(100);
  Vector2 尺寸 = {200, 200};

  Rectangle 单格1 = {0, 0, 尺寸.x, 尺寸.y};
  Rectangle 单格2 = {200, 0, 尺寸.x, 尺寸.y};
  Rectangle 单格3 = {0, 200, 尺寸.x, 尺寸.y};
  Rectangle 单格4 = {200, 200, 尺寸.x, 尺寸.y};

  // 主循环
  while (!WindowShouldClose()) // 检测窗口关闭按钮或ESC键
  {
    // Update
    //----------------------------------------------------------------------------------
    // TODO: 更新你的变量和逻辑
    //----------------------------------------------------------------------------------

    // 绘制
    //----------------------------------------------------------------------------------
    BeginDrawing();
    ClearBackground(GRAY);

    for (const auto& 点 : 泊松点集) {
      DrawCircleV({点.x * 尺寸.x, 点.y * 尺寸.y}, 2, BLACK);
    }
    DrawRectangleLinesEx(单格1, 3, BLACK);
    DrawText("1", 单格1.x + 5, 单格1.y + 5, 20, BLACK);

    for (const auto& 点 : 抖动网格点集) {
      DrawCircleV({点.x * 尺寸.x + 单格2.x, 点.y * 尺寸.y + 单格2.y}, 2, BLUE);
    }
    DrawRectangleLinesEx(单格2, 3, BLUE);
    DrawText("2", 单格2.x + 5, 单格2.y + 5, 20, BLUE);

    for (const auto& 点 : Vogel点集) {
      DrawCircleV({点.x * 尺寸.x + 单格3.x, 点.y * 尺寸.y + 单格3.y}, 2, RED);
    }
    DrawRectangleLinesEx(单格3, 3, RED);
    DrawText("3", 单格3.x + 5, 单格3.y + 5, 20, RED);

    for (const auto& 点 : Hammersley点集) {
      DrawCircleV({点.x * 尺寸.x + 单格4.x, 点.y * 尺寸.y + 单格4.y}, 2, GREEN);
    }
    DrawRectangleLinesEx(单格4, 3, GREEN);
    DrawText("4", 单格4.x + 5, 单格4.y + 5, 20, GREEN);

    EndDrawing();
    //----------------------------------------------------------------------------------
  }

  //--------------------------------------------------------------------------------------
  // 释放资源
  //
  //--------------------------------------------------------------------------------------
  CloseWindow(); // 关闭窗口和openGl上下文
  //--------------------------------------------------------------------------------------

  return 0;
}
