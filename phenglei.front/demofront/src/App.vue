<template>
  <div class="layout">
    <!-- 左侧菜单栏 -->
    <el-menu
        class="sidebar"
        mode="vertical"
        router
        default-active="/"
    >
      <el-menu-item index="/">参数配置</el-menu-item>
      <el-menu-item index="/build">构建运行</el-menu-item>
    </el-menu>

    <!-- 右侧内容 -->
    <div class="main-content">
      <router-view />
    </div>
  </div>
</template>

<style>
.layout {
  display: flex;
  height: 100vh;
}

.sidebar {
  width: 200px; /* 左侧菜单宽度 */
  min-width: 200px;
  border-right: 1px solid #ebeef5;
}

.main-content {
  flex: 1;
  padding: 16px;
  overflow-y: auto;
}
</style>