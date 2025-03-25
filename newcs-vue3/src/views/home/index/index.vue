<template>
<!-- 主页面右边 -->
<header class="header">
  <span class="title">欢迎您使用智能化学合成平台</span>
  <p class="description">
    本平台为基于强化学习的合成路线设计与平台，可以构建逆合成路线设计方法，形成人在回路下的知识引导路径优化方法，您也可以独立使用各功能模块完成独立的需求
  </p>
  <div class="underline"></div>
</header>

<div class="container">
  <div class="list-box">
    <router-link to="/XinxiChouqu" class="function-btn">
      <div class="card-content">
        <p class="card-title">反应参数信息抽取</p>
        <p class="card-subtitle">通过文献，化学方程等抽取化学反应</p>
        <div class="card-bottom">
          <img src="./image/xinxichouqu.png" alt="Icon1" class="card-image" />
          <router-link to="/help" class="card-link">帮助文档→</router-link>
        </div>
      </div>
    </router-link>
  </div>

  <div class="list-box">
    <router-link to="/ChanwuYuce" class="function-btn">
      <div class="card-content">
        <p class="card-title">反应产物结构预测</p>
        <p class="card-subtitle">给出反应物和反应条件，预测出产物结构</p>
        <div class="card-bottom">
          <img src="./image/chanwujiegou.png" alt="Icon2" class="card-image" />
          <router-link to="/help" class="card-link">帮助文档→</router-link>
        </div>
      </div>
    </router-link>
  </div>

  <div class="list-box">
    <router-link to="/LeixingYuce" class="function-btn">
      <div class="card-content">
        <p class="card-title">合成反应类型预测</p>
        <p class="card-subtitle">给出完整的化学反应，输出反应类型</p>
        <div class="card-bottom">
          <img src="./image/fanyingleixing.png" alt="Icon3" class="card-image"/>
          <router-link to="/help" class="card-link">帮助文档→</router-link>
        </div>
      </div>
    </router-link>
  </div>

  <div class="list-box">
    <router-link to="/ChanlvYuce" class="function-btn">
      <div class="card-content">
        <p class="card-title">化学反应产率预测</p>
        <p class="card-subtitle">给出完整的化学反应，输出产率</p>
        <div class="card-bottom">
          <img src="./image/fanyingchanlv.png" alt="Icon4" class="card-image"/>
          <router-link to="/help" class="card-link">帮助文档→</router-link>
        </div>
      </div>
    </router-link>
  </div>

  <div class="list-box">
    <router-link to="/TiaojianYouhua" class="function-btn">
      <div class="card-content">
        <p class="card-title">化学合成条件预测</p>
        <p class="card-subtitle">给出完整的化学反应，预测反应条件</p>
        <div class="card-bottom">
          <img src="./image/tiaojianyouhua.png" alt="Icon5" class="card-image" />
          <router-link to="/help" class="card-link">帮助文档→</router-link>
        </div>
      </div>
    </router-link>
  </div>

  <div class="list-box">
    <router-link to="/LujingYuce" class="function-btn">
      <div class="card-content">
        <p class="card-title">合成路径智能预测</p>
        <p class="card-subtitle">输入期望的产物，智能预测合成路径</p>
        <div class="card-bottom">
          <img src="./image/hechenglujing.png" alt="Icon6" class="card-image" />
          <router-link to="/help" class="card-link">帮助文档→</router-link>
        </div>
      </div>
    </router-link>
  </div>
</div>
  <!-- <div class="quick-start">
    <button class="quick-start-btn">查看帮助文档，快速上手</button>
  </div> -->
</template>
  

<script setup>
import { useRouter } from 'vue-router'
import { reactive } from 'vue'
const router = useRouter();  // 在组件内使用 router 对象来进行导航
</script>

<style scoped>
hr {
  border: none;             /* 移除默认边框 */
  height: 10px;             /* 设置高度 */
  background-color: #f70d0d;  /* 设置颜色 */
}

.header {
  height: 330px; /* 设置区域高度 */
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
}

.title {
  font-size: 50px; /* 标题字体大小 */
  color: #000; /* 标题字体颜色 */
  font-weight: bold; /* 标题加粗 */
}

.description {
  margin-top: 80px; /* 上方留出80px的空白 */
  font-size: 25px; /* 正文字体大小 */
  color: #000; /* 正文字体颜色 */
  font-weight: normal; /* 不加粗 */
  text-align: center;
  padding: 0 50px; /* 左右边距 */
}

.underline {
  width: calc(100% - 100px); /* 下划线的宽度为容器宽度的 100% 减去左右各 50 像素的边距 */
  height: 20px;
  background-color: rgb(224, 224, 224); /* 下划线颜色 */
  margin-top: 20px; /* 可根据需求调整距离 【underline的位置】距离【上面框】的距离为20px */
}

.container {
  display: flex;
  flex-wrap: wrap; /* 允许换行 */
  gap: 20px; /* 模块之间的间距 */
}
.list-box {
  /* width: 720px; */
  width: 530px; 
  flex: 0 0 auto; /* 保持固定宽度，不自动缩放 */
  margin-left: 60px;
  margin-bottom: 80px
}

.function-btn {
  display: block;
  text-decoration: none;
}

.card-content {
  position: relative; /* 相对定位，便于内部元素绝对定位 */
  padding: 50px;
  background-color: #fff;
  border: 1px solid #ddd;
  transition: background-color 0.3s ease;
  width: 100%;
  height: 100%;
  box-sizing: border-box;
}

.card-content:hover {
  background-color: rgb(229, 229, 229); /* 鼠标悬停时的背景颜色 */
}

.card-title {
  font-size: 35px;
  color:rgb(1,70,138);/* 标题颜色 */
  margin: 0;
}

.card-subtitle {
  font-size: 25px;
  font-weight: normal; /* 不加粗 */
  color: rgb(22, 22, 22); /* 副标题颜色 */
  margin: 20px 0;
}

.card-image {
  width: 110px;
  height: 110px;
  margin-top: 20px;
}

.card-link {
  position: absolute;
  bottom: 50px; /* 距离卡片底部50px */
  right: 50px; /* 距离卡片右侧50px */
  font-size: 25px;
  color: rgb(1,70,138); /* 链接颜色 */
  text-decoration: none;
}

.card-link:hover {
  text-decoration: underline; /* 鼠标悬停时的下划线 */
}


</style>
