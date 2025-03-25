<template>
  <div>
    <el-header class="header">
        合成路线智能预测
    </el-header>

    <el-tabs v-model="activeTab" class="tab-container" v-if = "showRoutes_db === 1">  <!-- 【标签页】组件，用于组织不同的内容面板 -->
      <el-tab-pane label="创建任务" name="create">
        <div class="container">
          <div class="column">
            <!-- 上传信息模块 -->
            <div class="section upload-section">  <!-- upload-section ? -->
              <h3>上传信息</h3>
              <div style="font-size: 1em;">示例分子:<br>
                  <p style="margin-left: 2em">O=C1C=CC2=C(CC[C@H](CC)[C@@H]2O)O1</p>
                  <p style="margin-left: 2em">Cc1cccc(c1N(CC(=O)Nc2ccc(cc2)c3ncon3)C(=O)C4CCS(=O)(=O)CC4)C</p>
              </div>
                <!-- <span>SMILES字符串:</span> -->
                <el-input v-model="moleculeInput" placeholder="输入SMILES字符串"  
                  clearable style="width: 52%;"></el-input><br>
                <!-- <span>期望路线:</span> -->
                <el-input v-model="number_of_routes" placeholder="输入期望的路线数量"  
                  clearable style="width: 52%;"></el-input><br>
                
              <el-button type="primary" class="submit-button" @click="submitTask">开始预测</el-button>
            </div>


            <!-- 预测结果模块 -->
            <div class="section prediction-results" v-if="predictionResults">
              <h3>预测结果</h3>
                  <hr>
                  全体图片个数:{{ totalImages_f }}<br>
                  路线最高分：{{ top_score_f }}<br>
                  <hr>
                  <div>
                    目标分子：
                    <img  :src="target_image_f" alt="目标分子结构图" class="target_image"/><br>
                  </div>

                  <hr>
                  <div v-for="(imgURL, index) in list" :key="index" class="reaction-item" >
                    <div class="reaction-info">
                      <h4>route:{{ index }}</h4>
                      <img @click="toggleExpand(index)" :src="imgURL" alt="分子结构图"  :class="['small_largeElement', { active: expandedIndex === index }]"/><br>
                    </div>
                    <el-button @click="toggleExpand(index)" class="expandButton">
                      {{ expandedIndex === index ? '折叠' : '展开' }}反应底物大图
                    </el-button>
                  </div>
                 <hr>
            </div>
          </div>
        </div>
      </el-tab-pane>

      <el-tab-pane label="最近结果" name="recent">
        <h1>最近结果</h1>
        <!-- 最近结果的内容，可以根据需求添加 -->



         <div class="recent-results" v-if="recentResults.length > 0">
          <div v-for="(result, index) in recentResults" :key="index" class="result-item">
            <h4>目标分子 ID: {{ index+1 }}</h4>
            <!-- result.id -->
            <img 
                :src="`data:image/png;base64,${result.b64_data}`" 
                alt="最近结果图像" 
                class="recent_tarimage"
                @click="() => { console.log('Clicked ID:', result.id); fetchPredictedRoutes(result.id); }"
             />
                
          </div>
         </div>
         <!-- 在前一个 v-if 条件为 false 时,展示下面的内容 -->
         <div v-else>  
          <p>没有最近结果可显示。</p>
         </div>



      </el-tab-pane>

      <el-tab-pane label="测试001" name="test">
        <h1>test001</h1>
        <!-- 最近结果的内容，可以根据需求添加 -->
      </el-tab-pane>
    </el-tabs>




    <div v-if = "showRoutes_db === 2">
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
  </div>
</template>

<script setup>
  import { ref,watch } from 'vue';
  import axios from 'axios';

  const moleculeInput = ref('');
  const number_of_routes = ref(1)  // 路线个数
  const totalImages_f = ref(0)  // 图像总数
  const top_score_f = ref(0)  // 最高分
  const target_image_f = ref('')

  const activeTab = ref("create")
  const predictionResults = ref(null)
  const showRoutes_db = ref(1)  // 展示数据库的页面
  let expandedIndex = ref(null)


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
    
    }
  }catch (error) {
    console.error('获取路线失败:', error);
  }
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
  if (newTab === 'recent') {
    fetch_targetimg_Results(); // 当切换到“最近结果”时获取数据
    
  }
});



  async function submitTask(){  // async异步要求
      console.log('1111') 
      
      result = await axios.post('http://127.0.0.1:3000/api/generate-molecule',
      {
        // 发送的数据
        molecule:{molecule: String(moleculeInput.value)},  
        numb_routes: number_of_routes.value,
      }
)
      console.log(result)

        // 接受的数据
      totalImages_f.value = result.data.total_images
      top_score_f.value = result.data.top_score
      target_image_f.value = `data:image/png;base64,${result.data.target_image}`; 


      // 图像结果处理 
      if (result.data.images && result.data.images.length > 0) {
        list.value = [];
        for (i=0;i<number_of_routes.value;i++){
          imgBase64 = `data:image/png;base64,${result.data.images[i]}`;
          list.value.push(imgBase64)
        }
      }
      predictionResults.value = 1
      // 可以执行
      // console.log('响应数据:', result.data)
      // imageUrl.value = `data:image/png;base64,${result.data.images[0]}`;
  }

  function toggleExpand(index){
     expandedIndex.value = (expandedIndex.value === index) ? null : index;
  }

  // function small_OR_large(){
  //   isActive.value = !isActive.value;
  // }







  
</script>

<style scoped>

.recent_tarimage {
  max-width: 100%;
  height: auto;
  cursor: pointer; /* 显示手型光标 */
  position: relative; /* 确保图像处于正确的层级 */
  z-index: 5; /* 提高图像的层级 */
}

/* .small_largeElement {
    background-color: red;
    width: 150px;
} */

/* .small_largeElement.active {
    background-color: blue;
    width: 200px;
} */



.header {
  background-color: rgb(1,70,138);
  padding: 0px;
  text-align: center;
  font-size: 1.5em;  /* 文字大小 */
  font-weight: bold;  /* 加粗 */
}
.tab-container {
  margin: 10px 0;
}
.container {
  display: flex;
  flex-wrap: wrap;
  gap: 10px;
  padding: 20px;
}
.column {
  /* flex: 1;
  允许 .column 元素在其父容器中占据可用的空间。
  在 Flexbox 布局中，这意味着该元素会根据其他兄弟元素的占用空间进行伸缩 */
  flex: 1;
  display: flex;
  flex-direction: column;
  
  gap: 10px;
}
.section {
  border: 1px solid #cccccc;
  padding: 10px;
  flex: 1;
  display: flex;              /* 使用 Flexbox */
  flex-direction: column;     /* 垂直排列子元素 */
  
}
.submit-button {
  margin-top: 10px;
  /* line-height: 40px;  */
  /* padding: 5px 10px;  调整内边距，减少按钮的高度和宽度 */
  font-size: 0.6em;    /* 调整字体大小 */
  /* min-width: 100px;    可选: 设置最小宽度 */
  max-width: 100px;
  height: 50px;
  
}
.reaction-item {  
  /* position: relative;                      使大图能够覆盖 */
  margin-bottom: 20px;
}
.reaction-info {
  
  display: flex;
  align-items: center;
}
.reaction-info h4 {
  margin-right: 10px;
}


.small_largeElement {
  margin-left: 1em;
  display: flex;
  flex-wrap: wrap;
  /* justify-content: center;    水平居中 */
  align-items: center; 
  width: 900px;
  /* height: auto; */
  /* height: 50px; */
}
.small_largeElement.active  {
  margin-left: 1em;
  width: 1800px;
  margin-top: 10px;
}

.recent_routesimage{
    margin-left: 1em;
    width: 1500px;
    margin-top: 10px;

}

.target_image{
  margin-left: 5em;
  width: 400px;
  max-width: 100%;        /* 确保图像不会超出容器 */
  display: flex;
  justify-content: center;  /* 水平居中 */
  align-items: center;      /* 垂直居中 */
  height: auto;             /* 根据需要调整高度 */
}

.expandButton{
  font-size: 0.7em;
  width: 200px;


}

</style>