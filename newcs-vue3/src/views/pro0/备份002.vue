<template>
  <el-container style="height: 100vh;">
    <!-- 顶部导航栏 -->
    <el-header class="header">
      测试0
    </el-header>
  <div>
    <h1>aiythfinder</h1>
    <h6>分子表达式:<input v-model="moleculeInput" placeholder="输入分子表达式" /></h6>

    <!-- type="number"代表只能输入数字     v-model.number代表字符5也可以变成数字5-->
    <h6>搜索路线个数:<input v-model.number="number_of_routes" placeholder="搜索路线个数" type="number" /></h6>

    <button @click="jumpTest">逆合成测试</button>
    <hr>
    
    <hr>
        O=C1C=CC2=C(CC[C@H](CC)[C@@H]2O)O1<br>
        O=C1C2([H])C=C(CO)[C@@H](O)[C@]3(O)[C@@H](O)C(C)=C[C@@]31[C@H](C)C[C@]4([H])[C@@]2([H])C4(C)C<br>
        Cc1cccc(c1N(CC(=O)Nc2ccc(cc2)c3ncon3)C(=O)C4CCS(=O)(=O)CC4)C<br>
    <hr>
    全体图片个数:{{ totalImages_f }}<br>
    路线最高分：{{ top_score_f }}<br>
    <div v-for="(imgURL, index) in list" :key="index">
        <img  :src="imgURL" alt="分子结构图" /><br>
        <span class="image-label">{{ index + 1 }}</span>
        <hr>
        <hr>
    </div>
    
    <hr>


  </div>

    
  </el-container>
</template>


<script setup>
  import { ref } from 'vue';
  import axios from 'axios';

  const moleculeInput = ref('');
  const number_of_routes = ref(1)  // 路线个数
  const totalImages_f = ref(0)  // 图像总数
  const top_score_f = ref(0)  // 最高分


  const imageUrl = ref('');
  let result = ref()
  let list = ref([])
  let i = 0
  let imgBase64 = 0

  async function jumpTest(){  // async异步要求
      console.log('1111') 
      
      result = await axios.post('http://127.0.0.1:3000/api/generate-molecule',
      {
        // 发送的数据
        molecule:{molecule: String(moleculeInput.value)},  
        numb_routes: number_of_routes.value,
      }
)
        // 接受的数据
      totalImages_f.value = result.data.total_images
      top_score_f.value = result.data.top_score

      // 图像结果处理 
      if (result.data.images && result.data.images.length > 0) {
        list.value = [];
        for (i=0;i<number_of_routes.value;i++){
          imgBase64 = `data:image/png;base64,${result.data.images[i]}`;
          list.value.push(imgBase64)
        }
      }

      // 可以执行
      // console.log('响应数据:', result.data)
      // imageUrl.value = `data:image/png;base64,${result.data.images[0]}`;
  }


</script>


<style scoped>
.header {
  background-color: rgb(1,70,138);
  color: white;
  text-align: center;
  padding: 0px 0;
}

.el-tabs {
  background-color: #ffffff;
  margin-bottom: 20px;
}

.module {
  padding: 20px;
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
  border-radius: 8px;
}

.el-upload__text {
  font-size: 12px;
  text-align: center;
}

.el-button {
  width: 100%;
}
</style>
