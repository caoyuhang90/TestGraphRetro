<template>
  <el-container class="layout-container-demo">
    <!-- 左侧菜单栏 -->
    <navMenu></navMenu>

    <el-container>
      <!-- 头信息 -->
      <!-- 动态面包屑，显示对应的功能名称 -->
      <el-header class="header-bar">
        <el-breadcrumb separator="/" class="crumbs">
          <el-breadcrumb-item v-for="(route, index) in $route.matched" :key="index" :to="{ path: route.path }">
            {{ route.meta.title }}
          </el-breadcrumb-item>
        </el-breadcrumb>
      </el-header>

      <!-- 内容区域 -->
      <el-main class="main-content">
        <RouterView />
      </el-main>
    </el-container>
  </el-container>
</template>

<script setup>
import navMenu from "./navMenu.vue";
</script>

<style scoped>
/* 主容器样式 */
.layout-container-demo {
  height: 100vh;
  display: flex;
  background-color: #f5f5f5;
}

/* 头部样式 */
.header-bar {
  position: relative; /* 相对于父元素定位 */
  background-color: var(--el-color-primary-light-7); /* 背景颜色 */
  color: var(--el-text-color-primary);  /* 文本颜色 */
}

.crumbs {
  padding: 10px 20px;
  font-size: 30px;
}

/* 菜单样式 */
.layout-container-demo .el-aside {
  background-color: #f0f2f5;
  color: var(--el-text-color-primary);
  border-right: 1px solid #e0e0e0;
  box-shadow: 2px 0 8px rgba(0, 0, 0, 0.1); /* 左侧菜单阴影 */
  padding-top: 20px;
}

.layout-container-demo .el-menu {
  border-right: none;
}

/* 主内容区域样式 */
.main-content {
  flex: 1;
  padding: 20px;
  background-color: #fff;
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1); /* 主内容区域阴影 */
  overflow-y: auto;
  border-radius: 8px;
  margin: 20px;
}

/* 工具栏样式 */
.toolbar {
  display: inline-flex;
  align-items: center;
  justify-content: flex-end;
  width: 100%;
}
</style>
