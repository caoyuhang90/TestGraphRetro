<template>
  <div class="common-layout">
    <el-container>
      <el-header class="header-bar">
        智能化学合成平台
      </el-header>
      <el-container class="layout-container-demo">
        <!-- 左侧 -->
        <navMenu></navMenu>  
        <el-main class="main-content">
          <!-- 右侧主内容区域,并包含了 <RouterView> 组件,用于渲染路由匹配的组件 -->
         <RouterView />
        </el-main>
      </el-container>
    </el-container>
  </div>
</template>

<script setup>
import navMenu from "./navMenu.vue";
</script>

<style scoped>
/* 主容器样式 */
.layout-container-demo {
  height: 100vh;
  display: flex;
}

.main-content{
  background-color: rgb(244, 244, 244)
}

/* 头部样式 */
.header-bar {
  position: relative; /* 相对于父元素定位 */
  background-color: rgb(60, 63, 65); /* 背景颜色 */
  color:white;  /* 文本颜色 */
  height: 70px;
  font-size: 24px;
  display: flex;
  align-items: center; /* 垂直居中 */
  text-indent: 2em; /* 缩进两格 */
}

/* 工具栏样式 */
.toolbar {
  display: inline-flex;
  align-items: center;
  justify-content: flex-end;
  width: 100%;
}
</style>
