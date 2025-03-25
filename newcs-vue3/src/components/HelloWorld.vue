<script setup>
import { ref } from 'vue'

// 定义一个属性，父组件可以通过 msg 属性向该组件传递字符串
defineProps({
  msg: String,
})

// 定义响应式变量 count
const count = ref(0)
</script>


<template>
<!-- 传入上面的msg属性 -->
  <h1>{{ msg }}</h1>

  <div class="card">
    <button type="button" @click="count++">count is {{ count }}</button>
    <p>
      Edit
      <code>components/HelloWorld.vue</code> to test HMR
    </p>
  </div>

  <!-- 链接到vue文档 -->
  <p>
    Check out
    <a href="https://vuejs.org/guide/quick-start.html#local" target="_blank"
      >create-vue</a
    >, the official Vue + Vite starter
  </p>
  <!-- 链接到IDE支持 -->
  <p>
    Learn more about IDE Support for Vue in the
    <a
      href="https://vuejs.org/guide/scaling-up/tooling.html#ide-support"
      target="_blank"
      >Vue Docs Scaling up Guide</a
    >.
  </p>
  <p class="read-the-docs">Click on the Vite and Vue logos to learn more</p>
</template>

<!-- 定义一些样式 -->
<style scoped>
.read-the-docs {
  color: #888;
}
</style>
