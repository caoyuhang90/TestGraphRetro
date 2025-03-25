<template>
  <el-aside class="menu-aside" width="400px">
    <el-scrollbar>  <!-- 滚动条 -->    <!-- -->
      <!-- <el-menu>: 用于创建一个可折叠的菜单。
           default-openeds 属性定义【初始展开】的子菜单项 . 要不要展开
           router 属性表示使用 Vue Router 进行导航。 -->
      <el-menu :default-openeds="['0', '1', '2','3']" router class="menu-list">
        <!-- 遍历 items 数组生成菜单 -->

        <!-- 大菜单 -->
        <el-sub-menu  
          v-for="(item, index) in items"
          :key="index"
          :index="index.toString()"  
          class="sub-menu"
        >  <!-- 以确保index在菜单中唯一  -->
        
          <template #title>    
            <el-icon class="menu-icon"><Menu /></el-icon>  <!-- 图标 -->
            <span class="menu-text">{{ item.text }}--{{ index.toString() }}</span>  <!-- span 行内标签 -->
          </template>

          <!-- 遍历子菜单项，将 index 绑定为 url,用作【路由导航】 -->
          <el-menu-item
            v-for="(child, idx) in item.child"
            :key="idx"
            :index="child.url"  
            class="menu-item"  
          > 
            <el-icon class="submenu-icon"><Menu /></el-icon>
            {{ child.text }}--{{ child.url }}
          </el-menu-item>
        </el-sub-menu>
      </el-menu>
    </el-scrollbar>
  </el-aside>
</template>

<script setup>
import { reactive } from 'vue'
import { Menu } from '@element-plus/icons-vue'  // 图标

const items = reactive([
  {
    text: '测试001',
    url: '',
    child: [
      { text: 'test', url: '/test0' },
      // { text: '项目管理', url: '/project' },
      // { text: '多人合作', url: '/Collaborators' },
    ],

  },


  {
    text: '工作区',
    url: '',
    child: [
      { text: '个人主页', url: '/HomePage' },
      { text: '项目管理', url: '/project' },
      { text: '多人合作', url: '/Collaborators' },
    ],
  },

  {
    text: '功能区',
    url: '',
    child: [
      { text: '信息抽取', url: '/XinxiChouqu' },
      { text: '产物预测', url: '/ChanwuYuce' },
      { text: '条件预测', url: '/TiaojianYouhua' },
      { text: '类型预测', url: '/LeixingYuce' },
      { text: '产率预测', url: '/ChanlvYuce' },
      { text: '路径规划', url: '/LujingYuce' },
      { text: '自选断键', url: '/ZiXuanDuanJian' },
    ],
  },

  {
    text: '更多信息',
    url: '',
    child: [
      { text: '快速上手', url: '/quick' },
      { text: '关于我们', url: '/about' },
    ],
  },
])

</script>

<style scoped>
.menu-aside {
  background-color: #f5f5f5;
  height: 100vh;
  color: #000;
  border-right: 1px solid #e0e0e0;
}

.menu-text {
  color: #000;
  font-weight: 600;
  padding-left: 10px;
  font-size: 24px;
}

.menu-item {
  display: flex;
  align-items: center;
  color: #000;
  font-size: 24px;
  padding-left: 20px;
  height: 70px;
}

.el-menu-item.is-active {
  background-color: rgb(228,228,228);
  color: rgb(22, 22, 22);
}

.el-scrollbar__wrap {
  scrollbar-width: thin;
  scrollbar-color: #2c3e50 #f0f0f0;
}

.el-scrollbar__wrap::-webkit-scrollbar {
  width: 8px;
}

.el-scrollbar__wrap::-webkit-scrollbar-track {
  background-color: #2c3e50;
}

.el-scrollbar__wrap::-webkit-scrollbar-thumb {
  background-color: #3498db;
  border-radius: 4px;
}

.menu-list {
  width: 400px;
  height: 100vh;
}

.menu-icon,
.submenu-icon {
  margin-right: 40px;
  font-size: 24px;
}
</style>
  
  
  
  
<!-- <template>
<el-aside width="250px">
        
    <el-scrollbar >
      <el-menu :default-openeds="['1', '3']">

        <el-sub-menu index="1">
          <template #title>
            <el-icon><message /></el-icon>工作区
          </template>
            <el-menu-item index="1-1">个人主页</el-menu-item>
            <el-menu-item index="1-2">项目管理</el-menu-item>
            <el-menu-item index="1-3">多人合作</el-menu-item>          
        </el-sub-menu>

        <el-sub-menu index="2">
          <template #title>
            <el-icon><icon-menu /></el-icon>功能区
          </template>
            <el-menu-item index="2-1">信息抽取</el-menu-item>
            <el-menu-item index="2-2">产物预测</el-menu-item>
            <el-menu-item index="2-3">条件优化</el-menu-item>
            <el-menu-item index="2-4">类型预测</el-menu-item>
            <el-menu-item index="2-5">产率预测</el-menu-item>
            <el-menu-item index="2-6">路径规划</el-menu-item>
        </el-sub-menu>

        <el-sub-menu index="3">
          <template #title>
            <el-icon><setting /></el-icon>更多信息
          </template>
            <el-menu-item index="3-1">快速上手</el-menu-item>
            <el-menu-item index="3-2">关于我们</el-menu-item>
        </el-sub-menu>
      </el-menu>
    </el-scrollbar>
  </el-aside>
</template>

<script setup>
import {useRouter,useRoute} from 'vue-router'
import { reactive } from 'vue'

const items = reactive([
    {
        text: '工作区',
        url:'',
        child:[
        { text:'个人主页', url:'/newhome'},
        { text:'项目管理', url:'/project'},
        { text:'多人合作', url:'/Collaborators'}
    ]
    },

    {
        text: '功能区',
        url:'',
        child:[
        { text:'信息抽取', url:'/XinxiChouqu'},
        { text:'产物预测', url:'/ChanwuYuce'},
        { text:'条件优化', url:'/TiaojianYouhua'},
        { text:'类型预测', url:'/LeixingYuce'},
        { text:'产率预测', url:'/ChanlvYuce'},
        { text:'路径规划', url:'/LujingYuce'}
    ]
    },

    {
        text: '更多信息',
        url:'',
        child:[
        { text:'快速上手', url:'/quick'},
        { text:'关于我们', url:'/about'}
    ] 
    }
])
</script> -->