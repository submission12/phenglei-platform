<template>
  <div class="container">

    <h2 class="title">计算算例</h2>

    <!-- ================= 构建 ================= -->
    <div class="block">

      <el-button @click="buildEnvExe" :loading="buildRunning">
        构建当前环境可执行文件
      </el-button>

      <el-progress
          v-if="buildRunning || buildSuccess"
          :percentage="buildProgress"
          :status="buildSuccess ? 'success' : undefined"
          style="margin-top: 10px;"
      />

    </div>

    <!-- ================= 运行 ================= -->
    <div class="block">

      <el-button
          @click="runCase"
          :loading="runRunning"
          :disabled="runRunning"
      >
        运行算例（实时）
      </el-button>

      <div class="hint">
        <span v-if="runRunning">正在运行算例（实时监控中）...</span>
        <span v-else-if="runSuccess">✔ 运行完成</span>
      </div>

    </div>

    <!-- ================= 实时日志 ================= -->
    <div class="block">

      <div class="log-box">
        <div v-for="(line, index) in logs" :key="index">
          {{ line }}
        </div>
      </div>

    </div>

  </div>
</template>

<script setup>
import { ref, onUnmounted } from "vue";
import axios from "axios";
import { ElMessage } from "element-plus";

/* ================= 构建状态 ================= */
const buildRunning = ref(false);
const buildProgress = ref(0);
const buildSuccess = ref(false);

/* ================= 运行状态 ================= */
const runRunning = ref(false);
const runSuccess = ref(false);

/* ================= 日志 ================= */
const logs = ref([]);

/* ================= SSE ================= */
let runES = null;

/* 🚨 关键：是否正常结束 */
let isFinished = false;

/* ================= SSE连接 ================= */
function connectRunSSE(taskId) {

  const es = new EventSource(
      `http://localhost:8080/log/subscribe/${taskId}`
  );

  // ===== 日志 =====
  es.addEventListener("log", (e) => {
    logs.value.push(e.data);

    setTimeout(() => {
      const box = document.querySelector(".log-box");
      if (box) box.scrollTop = box.scrollHeight;
    }, 0);
  });

  // ===== 完成 =====
  es.addEventListener("done", () => {

    isFinished = true;

    runRunning.value = false;
    runSuccess.value = true;

    ElMessage.success("运行完成");
  });

  // ===== 错误 =====
  es.onerror = () => {

    runRunning.value = false;
    es.close();

    // 🚨 关键：如果已经正常结束，就不报错
    if (isFinished) return;

    ElMessage.error("SSE连接异常");
  };

  return es;
}

/* ================= 运行算例 ================= */
async function runCase() {

  const taskId = "run";

  // 防止重复连接
  if (runES) {
    runES.close();
    runES = null;
  }

  runRunning.value = true;
  runSuccess.value = false;
  logs.value = [];

  // 🚨 每次运行重置状态
  isFinished = false;

  runES = connectRunSSE(taskId);

  try {
    await axios.get("http://localhost:8080/run-exe");
    ElMessage.info("算例已启动，实时监控中...");
  } catch (e) {
    runRunning.value = false;
    runES?.close();
    ElMessage.error("运行失败");
  }
}

/* ================= 构建 ================= */
async function buildEnvExe() {
  try {
    buildRunning.value = true;

    await axios.get("http://localhost:8080/build-exe");

    buildProgress.value = 100;
    buildSuccess.value = true;
    buildRunning.value = false;

    ElMessage.success("构建成功");

  } catch (e) {
    buildRunning.value = false;
    ElMessage.error("构建失败");
  }
}

/* ================= 清理 ================= */
onUnmounted(() => {
  runES?.close();
});
</script>

<style scoped>
.container {
  padding: 20px;
}

.title {
  font-size: 20px;
  font-weight: 600;
  margin-bottom: 20px;
}

.block {
  margin-bottom: 30px;
}

.hint {
  margin-top: 6px;
  color: #666;
  font-size: 13px;
}

.log-box {
  height: 250px;
  overflow-y: auto;
  background: #000;
  color: #00ff00;
  padding: 10px;
  font-family: monospace;
  font-size: 12px;
  border-radius: 6px;
}
</style>