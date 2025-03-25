import { createApp } from 'vue'
import './style.css'

import ElementPlus from 'element-plus'  //导入ElementPlus所有模块和功能
import 'element-plus/dist/index.css'  // 导入ElementPlus所需的css

import * as ElementPlusIconsVue from '@element-plus/icons-vue'


import App from './App.vue'

import router from './router'



// createApp(App)
//   .use(router)      // 先注册 router
//   .use(ElementPlus) // 再注册 ElementPlus
//   .mount('#app');   // 最后挂载应用



const app = createApp(App)
app.use(router)
  .use(ElementPlus)

for (const [key, component] of Object.entries(ElementPlusIconsVue)) {
  app.component(key, component)
}

app.mount('#app')