#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

class Etcher
{
public:
    // 构造：设置默认位姿、长度与最大倾角（弧度）
    Etcher();

    // 位姿与参数（公开以便 UI 直接读写）
    glm::vec3 position; // world-space position (x,y,z)
    glm::vec3 rotation; // rotation angles (rx, ry, rz) in radians
    float length;       // length along local -Z (used for centering)
    float maxTiltX;
    float maxTiltY;
    float maxTiltZ;

    // 限制 rotation 在允许范围内（在 UI 回调后调用）
    void clampRotation();

    // 计算并返回刻蚀头在世界空间的模型矩阵（与 Viewer::updateEtcherMesh 中逻辑一致）
    glm::mat4 transform() const;

    // 便捷设置
    void setPosition(const glm::vec3 &p) { position = p; }
    void setRotation(const glm::vec3 &r)
    {
        rotation = r;
        clampRotation();
    }
};

// 在单独的翻译单元中定义一个全局实例，便于现有代码直接使用 `etcher`
extern Etcher etcher;