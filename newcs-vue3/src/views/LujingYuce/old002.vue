<template>
  <!-- style="height: 100vh;" -->
  <el-container >
    <el-header class="header">
      合成路线智能预测
    </el-header>

    <p style="margin-left: 2em">O=C1C=CC2=C(CC[C@H](CC)[C@@H]2O)O1</p>

    <!-- flex-direction: column:设置 Flexbox 的主轴方向为纵向。这意味着子元素将垂直排列 -->
    <el-main style="display: flex; flex-direction: column; height: 100%;">  
      <!-- <p v-for="o in 4" :key="o" class="text item">{{ 'List item ' + o }}</p> -->
      <!-- 卡片一 -->
      <el-card class="upload-card">
        <h3 class="title">上传信息</h3>
        <!-- 分栏布局 -->
         <div >
          <el-row :gutter="20" class="button-row">
            <el-button type="primary" @click="activeTab = 'smiles'">smiles反应式</el-button>
            <el-button type="primary" @click="activeTab = 'histryRoutes'">历史记录</el-button>
          </el-row>
          <!-- <el-row :gutter="20" class="button-row">
            <el-button type="primary" @click="activeTab = 'histryRoutes'">历史记录</el-button>
          </el-row> -->
         </div>

        <!-- 表单 
         label-width：表单项标签的宽度为 100 像素
         @submit.prevent：当表单提交时，调用 getPredictedYield 方法，并阻止默认的提交行为
         -->  
        <el-form v-if="activeTab === 'smiles'"  :model="form" label-width="150px"  @submit.prevent="submitTask">
          <el-form-item label="smiles">
            <el-input v-model="form.smiles" placeholder="输入化学反应的 SMILES 表达式"></el-input>
          </el-form-item>
          <el-form-item label="期望路线个数">
            <el-input type='number' v-model.number="form.number_of_routes" placeholder="输入期望的路线"></el-input>
          </el-form-item>

          <el-form-item>
            <el-button type="primary" @click="()=>{submitTask()}">提交</el-button>
            <el-icon  v-if="loading_flag === 1" class="is-loading" style="margin-left: 2em;font-size: 25px;" ><Loading /></el-icon>  <!-- loading图标 -->
          </el-form-item>
        </el-form>
      </el-card>


      
      <el-card v-if="activeTab === 'smiles'"  class="result-card  " :style="{ flex: predictionResults === 0 ? '1' : 'unset' }">
        <h2 class="title" style="font-size: 32px !important ">合成路线</h2>

        <!-- <p v-for="o in 4" :key="o" class="text item">{{ 'List item ' + o }}</p> -->
        <div v-if="predictionResults">
            <el-row :gutter="10">
              <el-col :span="2" >
                <span class="txt_target">目标分子：</span>
              </el-col>
              
              <el-col :span="22">
                <img  :src="target_image_f" alt="目标分子结构图" class="target_image"/>
              </el-col>
            </el-row>



            <el-row :gutter="10">
              <el-col :span="2" >
                <span style="font-size: 20px; width: 100px; ">选择路线:</span>
              </el-col>
              
              <el-col :span="22">
                <el-select  class="slectRoutes" v-model="selectedImage" placeholder="选择结果图"   style="width: 300px;"  >
                  <el-option v-for="(score, index) in all_scores_f" :key="index" 
                  :label="'路线 ' + (index + 1) + '&nbsp;&nbsp&nbsp;&nbsp&nbsp;&nbsp&nbsp;&nbsp&nbsp;&nbsp得分' + score" :value="index"></el-option>
                </el-select>
              </el-col>
            </el-row>

          <div v-if="predictionResults">
            <!-- 显示图像 -->
            <img :src="list[selectedImage]" alt="result image" style="max-width: 100%; height: auto; padding: 10px;">
          </div>
        </div>
        <p v-else class="result-title">请输入 产物SMILES 表达式并提交以查看预测结果。</p>
      </el-card>



      <!-- 历史记录 -->
      <el-card v-if="activeTab === 'histryRoutes' && showRoutes_db === 1" class="result-card  ">
         <div class="recent-results" v-if="recentResults.length > 0">
          <div v-for="(result, index) in recentResults" :key="index" class="result-item">



              <el-card class="result-card-target">
                <p style="margin: 10px;font-size: 20px;">目标分子 ID: {{ index+1 }}</p>
                <p style="margin: 10px;font-size: 20px;;width: 300px;" >创建时间: {{ formatDate(result.time) }}</p>
                <!-- <p>日期: <el-date-picker v-model="result.time" type="datetime" format="yyyy-MM-dd HH:mm:ss" /></p> -->

                <!-- result.id -->
                <el-image 
                    :src="`data:image/png;base64,${result.b64_data}`" 
                    alt="最近结果图像" 
                    class="recent_tarimage"
                    @click="() => { console.log('Clicked ID:', result.id); fetchPredictedRoutes(result.id); }"
                />
              <el-icon @click="delet_targetimg(result.id)" style="cursor: pointer;"><Delete /></el-icon>
              <!-- <el-button type="primary" @click="()=>{delet_targetimg(result.id),fetch_targetimg_Results()}">删除</el-button> -->

              </el-card>

          </div>
         </div>
         <!-- 在前一个 v-if 条件为 false 时,展示下面的内容 -->
         <div v-else>  
          <p>没有最近结果可显示。</p>
         </div>
        


      </el-card >
          

      <el-card v-if = "activeTab === 'histryRoutes' && showRoutes_db === 2"  class="result-card">

        <div >
          <h1>数据库路线展示</h1>
          <el-button type="primary" class="submit-button" @click="Back_tabs">返回</el-button>

          <div v-for="(routes, index) in predictRoutes" :key="index" class="result-item">
            <h4>路线 ID: {{ index+1 }}</h4>
            <!-- routes.id -->
            <img 
                :src="`data:image/png;base64,${routes.b64_routesdata}`" 
                alt="路线" 
                class="recent_routesimage"
              />
          </div>
        </div>

      </el-card>

    </el-main>
  </el-container>
</template>

<script setup >
  import { ref,watch } from 'vue';
  import axios from 'axios';

  const moleculeInput = ref('');

  let loading_flag = ref(0) // lodaing符号得使用标志


  // const number_of_routes = ref(1)  // 路线个数
  const totalImages_f = ref(0)  // 图像总数
  const top_score_f = ref(0)  // 最高分
  const target_image_f = ref('')
  const all_scores_f = ref([])  // 所有的分数

  const activeTab = ref("smiles")
  const predictionResults = ref(null)
  const showRoutes_db = ref(1)  // 展示数据库的页面
  let expandedIndex = ref(null)

  // const activeTab= ref('smiles')



  const form = ref({
    smiles: '' ,// 表单：初始化 smiles 字段
    number_of_routes: null
  })

  const selectedImage = ref(0)





  let result = ref()
  let list = ref([])
  const recentResults = ref([]);  
  const predictRoutes = ref([]);
  let i = 0
  let imgBase64 = 0

function Back_tabs(){

      predictRoutes.value = []    // 重置
      showRoutes_db.value = 1

}




async function fetchPredictedRoutes(tarimg_id) {
  try{
    console.log('id:tarimg_id')
    const response = await axios.get(`http://127.0.0.1:3000/user/showRoutes/${tarimg_id}`);
    console.log('response',response)
    if (response.data) {
      console.log('routeIN')
      
      predictRoutes.value = response.data.data; // 保存最近结果
      console.log('routes',predictRoutes.value)
      showRoutes_db.value = 2
      console.log('showRoutes_db.value',showRoutes_db.value)
    }
  }catch (error) {
    console.error('获取路线失败:', error);
  }
}

// 删除目标图像   http://127.0.0.1:3000/user/delete-targetimg/37   http://127.0.0.1:3000/user/delete-targetimg/31
async function delet_targetimg(tarimg_id) {
  try{
    console.log('执行删除')
    const response = await axios.delete(`http://127.0.0.1:3000/user/delete-targetimg/${tarimg_id}`)


    // 删除后的更新操作
    recentResults.value = recentResults.value.filter(result => result.id !== tarimg_id);
  }catch(error){
    console.error('删除失败',error)
  }
  
}

function formatDate(dateString) {
  const date = new Date(dateString);
  const year = date.getFullYear();
  const month = String(date.getMonth() + 1).padStart(2, '0'); // 月份是从0开始的，所以要+1
  const day = String(date.getDate()).padStart(2, '0');
  const hours = String(date.getHours()).padStart(2, '0');
  const minutes = String(date.getMinutes()).padStart(2, '0');
  const seconds = String(date.getSeconds()).padStart(2, '0');
  return `${year}-${month}-${day} ${hours}:${minutes}:${seconds}`;
}


  // 获取最近target图像 
async function fetch_targetimg_Results() {
  try {
    const response = await axios.get('http://127.0.0.1:3000/user/1/showTarget/');
    if (response.data && response.data.data) {
      recentResults.value = response.data.data; // 保存最近结果
      console.log(response)
      

      // response结构 {
      //   data{
              //     data{
              //       {id=xxx,b64_data=xxx}
              //       {id=xxx,b64_data=xxx}
              //       {id=xxx,b64_data=xxx}
              //     }
              //     message：
      //   }
      // }


    }
  } catch (error) {
    console.error('获取最近结果失败:', error);
  }
}

// 监听选项卡变化
watch(activeTab, (newTab) => {
  if (newTab === 'histryRoutes') {
    fetch_targetimg_Results(); // 当切换到“最近结果”时获取数据
    console.log('监听选项卡变化')
  }

  if (newTab === 'smiles') {
    showRoutes_db.value = 1
  }
});




  async function submitTask(){  // async异步要求

      
      console.log('loading_flag',loading_flag) 

      console.log(form)
      console.log(form.value.smiles)
      console.log(form.value.number_of_routes)
      loading_flag.value = 1
      result = await axios.post('http://127.0.0.1:3000/api/generate-molecule',
      {
        // 发送的数据
        molecule:{molecule: String(form.value.smiles)},  
        numb_routes: form.value.number_of_routes,
      }
)
      console.log(result)

        // 接受的数据
      totalImages_f.value = result.data.total_images  // 路线图像
      top_score_f.value = result.data.top_score  // 最高分数
      target_image_f.value = `data:image/png;base64,${result.data.target_image}`;   // 目标分子图像
      all_scores_f.value = result.data.all_scores  // 分数
      console.log('all_scores_f:',all_scores_f)
      // console.log('totalImages_f',totalImages_f)

      // 图像结果处理 
      if (result.data.images && result.data.images.length > 0) {
        list.value = [];
        for (i=0;i<form.value.number_of_routes;i++){
          imgBase64 = `data:image/png;base64,${result.data.images[i]}`;
          list.value.push(imgBase64)
          
        }
        console.log('list0',list)
        loading_flag.value = 0  // loading 图标出现
      }
      predictionResults.value = 1
      console.log('list1:',list)
      
      // 可以执行
      // console.log('响应数据:', result.data)
      // imageUrl.value = `data:image/png;base64,${result.data.images[0]}`;
  }

  function toggleExpand(index){
     expandedIndex.value = (expandedIndex.value === index) ? null : index;
  }

</script>

<style scoped>



.target_image{
  display: flex;               /* 使用 Flexbox */
  align-items: center;        /* 垂直居中对齐 */
  justify-content: center; 
  width: 300px;     

}


.txt_target {

  display: flex;
  justify-content: center; /* 水平居中 */
  align-items: center;    /* 垂直居中 */
  width: 100px;
  height: 300px;
  font-size: 20px;
}


.title{
  text-align: center;
}


.header {
  background-color:rgb(1,70,138);
  color: white;
  text-align: center;
  padding: 0px 0;
}

.upload-card {
  margin-top: 20px;
  padding: 20px;
  border: 1px solid #dcdfe6;
  border-radius: 8px;
}

.button-row {
  margin-bottom: 20px;
}


/* 卡片样式 */
.result-card-target {
  display: flex;               /* 使用 Flexbox */
  align-items: center;        /* 垂直居中对齐 */


  margin-top: 10px;
  padding: 25px;
  /* text-align: center; */
  border: 1px solid #dcdfe6;
  border-radius: 8px;  /* 圆角边框，使卡片看起来更柔和 */
  flex: 1; 
  width: 350px;
  height: 450px;
}

/* 鼠标悬停时的背景颜色 */
.result-card-target:hover {
  background-color: rgb(229, 229, 229); 
} 


/* 每个卡片 */
.result-card {     
  margin-top: 20px;
  padding: 20px;
  /* text-align: center; */
  border: 1px solid #dcdfe6;
  border-radius: 8px;  /* 圆角边框，使卡片看起来更柔和 */
  flex: 1; 
}


/* 鼠标悬停时的背景颜色 */
/* :hover {
  background-color: rgb(229, 229, 229); 
} */

.result-card p {
  font-size: 24px;
  color: #5D9CEC;
  font-weight: bold;
}
.result-title {
  text-align: center;

  font-size: 32px !important;
  color:rgb(1, 70, 138) !important;
}





:deep(.el-button--primary) {   /* 设置主按钮格式 */
  background-color: rgb(1, 70, 138) !important;
  color: #ffffff !important;
  border-color: rgb(1, 70, 138) !important;
}
:deep(.el-form-item__label) {  /* 设置表单字体 */
  font-size: 20px;
}

/* 调整 el-select 和 el-option 的字体大小 */
:deep(.el-select .el-input__inner),
:deep(.el-option) {
  font-size: 500px ;
}

/* 调整 el-option 的下拉菜单字体大小 */
:deep(.el-select-dropdown__item) {
  font-size: 500px ;
}


.recent_tarimage {
  max-width: 100%;
  
  height: auto;
  cursor: pointer; /* 显示手型光标 */
  position: relative; /* 确保图像处于正确的层级 */
  z-index: 5; /* 提高图像的层级 */
 margin-bottom: 20px;
}


.recent_routesimage{
  /* width: 95%; */
  max-width: 95%;
}


.recent-results{
  position: relative;
    display: flex;
  flex-wrap: wrap; /* 允许换行 */
  gap: 20px; /* 模块之间的间距 */
}


</style>
