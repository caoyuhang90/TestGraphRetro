<template>
  <div>
    <el-header class="header">
        自选断键
    </el-header>

    <div class="container">
      <div class="column">
        <div class="section ">  
              <h3>上传信息</h3>
              <p>C[C@@]12[C@@](C(SCF)=O)([C@@H](C[C@]1([C@@]3(C[C@@H](C4=CC(C=C[C@@]4([C@]3([C@H](C2)O)F)C)=O)F)[H])[H])C)OC(CC)=O</p>
              <p>15 16 28 30</p>
              <p>CCNC(=O)CCC/C=C\C[C@H]1[C@H](C[C@H]([C@@H]1/C=C/[C@H](CCC2=CC=CC=C2)O)O)O</p>
              <p>9 10 14 15</p>
              
              
              <!-- <div style="font-size: 1em;">示例分子:<br>
                  <p style="margin-left: 2em">O=C1C=CC2=C(CC[C@H](CC)[C@@H]2O)O1</p>
                  <p style="margin-left: 2em">Cc1cccc(c1N(CC(=O)Nc2ccc(cc2)c3ncon3)C(=O)C4CCS(=O)(=O)CC4)C</p>
              </div> -->

                <!-- <span>SMILES字符串:</span> -->
                <!-- 原子映射 -->
                
                 <el-card style="width: 90%;">
                  <div >
                    <div class="card1-layout1">
                      <el-input class="input_smiles" v-model="moleculeInput" placeholder="输入SMILES字符串"  
                      clearable style="width: 60%;"></el-input>
                      <el-button type="primary" class="submit-button" @click="TragetImage">原子映射</el-button>

                    </div>

                    <div>
                      <img 
                      v-if = 'target_image_f'
                      :src="`data:image/png;base64,${target_image_f}`" 
                      alt="最近结果图像" 
                      class="recent_tarimage"
                       />
                    </div>
                  </div>


                 </el-card>



                <!-- 合成子预览   :placeholder="`输入期望的第 ${index + 1} 组断键`"-->
                 <el-card style="width: 90%;">

                  <el-tabs v-model="activeTab">
                    <el-tab-pane label="单原子" name="danyuanzi"
                    style="display: flex; flex-direction: column; gap: 10px;">
                      <el-input class="input_smiles" v-model="oneAtom" placeholder="输入期望的单个原子序号"  
                    clearable style="width: 52%;"></el-input>
                     <el-input class="input_smiles" v-model="beam_num0" placeholder="束宽"  
                    clearable style="width: 20%;"></el-input>
                    <el-button type="primary" class="submit-button" @click="HeChengZi_oneAtmo">合成子预览</el-button>

              

                    </el-tab-pane>

                    <el-tab-pane label="单键" name="danjian" 
                    style="display: flex; flex-direction: column; gap: 10px;">
                      <el-input class="input_smiles" v-model="oneBond1" placeholder="输入期望的第一个原子"  
                    clearable style="width: 52%;"></el-input>
                      <el-input class="input_smiles" v-model="oneBond2" placeholder="输入期望的第二个原子"  
                    clearable style="width: 52%;"></el-input>
                    <el-input class="input_smiles" v-model="beam_num1" placeholder="束宽"  
                    clearable style="width: 20%;"></el-input>
                    <el-button type="primary" class="submit-button" @click="HeChengZi_oneBond">合成子预览</el-button>
                    </el-tab-pane>

                    <el-tab-pane label="双键" name="shuangjian"
                    style="display: flex; flex-direction: column; gap: 10px;">
                      <el-input class="input_smiles" v-model="first_edit1" placeholder="输入期望的第一组断键的第一个原子"  
                      clearable style="width: 52%;"></el-input>
                      <el-input class="input_smiles" v-model="first_edit2" placeholder="输入期望的第一组断键的第二个原子"  
                      clearable style="width: 52%;"></el-input>
                      <el-input class="input_smiles" v-model="second_edit1" placeholder="输入期望的第二组断键的第一个原子"  
                      clearable style="width: 52%;"></el-input>
                      <el-input class="input_smiles" v-model="second_edit2" placeholder="输入期望的第二组断键的第二个原子"  
                      clearable style="width: 52%;"></el-input>
                      <el-input class="input_smiles" v-model="beam_num2" placeholder="束宽1"  
                    clearable style="width: 20%;"></el-input>
                      <el-input class="input_smiles" v-model="beam_num22" placeholder="束宽2"  
                    clearable style="width: 20%;"></el-input>
                      <el-button type="primary" class="submit-button" @click="HeChengZi">合成子预览</el-button>
                    </el-tab-pane>

                  </el-tabs>








                  <!-- 增删
                  <div v-for="(input, index) in inputs" :key="index" class="input-item">
                    <el-input
                      v-model="input.value"
                      style="width: 10%;"
                    ></el-input>
                    <el-button
                      type="danger"
                      @click="removeInput(index)"
                      v-if="inputs.length > 1"
                    >
                      删除
                    </el-button>
                  </div>
                 <el-button type="primary" @click="addInput">添加原子序号</el-button> -->



                  <!-- <div class="card1-layout1">
                    <el-input class="input_smiles" v-model="first_edit" placeholder="输入期望的第一组断键"  
                    clearable style="width: 52%;"></el-input>
                    <el-input class="" v-model="second_edit" placeholder="输入期望的第二组断键"  
                    clearable style="width: 52%;"></el-input>
                    <el-button type="primary" class="submit-button" @click="HeChengZi">合成子预览</el-button>
                  </div> -->


                  <img 
                  v-if = "heCzi_img_danyuanzi &&  activeTab === 'danyuanzi'"
                  :src="`data:image/png;base64,${heCzi_img_danyuanzi}`" 
                  alt="最近结果图像" 
                  class="heCzi_image"
                 />

                  <img 
                  v-if = "heCzi_img_danjian &&  activeTab === 'danjian'"
                  :src="`data:image/png;base64,${heCzi_img_danjian}`" 
                  alt="最近结果图像" 
                  class="heCzi_image"
                 />

                  <img 
                  v-if = "heCzi_img &&  activeTab === 'shuangjian'"
                  :src="`data:image/png;base64,${heCzi_img}`" 
                  alt="最近结果图像" 
                  class="heCzi_image"
                 />

                 </el-card>


                <el-card style="width: 90%;">
                  <!-- 任务提交 -->
                  <el-button type="primary" class="submit-button" @click="edit_context">提交任务</el-button>
                  <!-- <img 
                  v-if="result_img "
                  :src="`data:image/png;base64,${result_img}`" 
                  alt="最近结果图像" 
                  class="recent_tarimage"
                  /> -->
                  <div v-if="recentResults0 && activeTab === 'danyuanzi'">
                    <div v-for=" (imgB64, index) in recentResults0" class="reaction-info">
                        <h4>结果:{{ index }} -- 得分：{{ scores_list0[index] }}</h4>
                        <img  :src="`data:image/png;base64,${imgB64}`" alt="分子结构图"  /><br>
                    </div>
                  </div>

                  <div v-if="recentResults1 && activeTab === 'danjian'">
                    <div v-for=" (imgB64, index) in recentResults1" class="reaction-info">
                        <h4>结果:{{ index }} -- 得分：{{ scores_list1[index] }}</h4>
                        <img  :src="`data:image/png;base64,${imgB64}`" alt="分子结构图"  /><br>
                    </div>
                  </div>


                  <div v-if="recentResults2 && activeTab === 'shuangjian'">
                    <div v-for=" (imgB64, index) in recentResults2" class="reaction-info">
                        <h4>结果:{{ index }} -- 得分：{{ scores_list2[index] }}</h4>
                        <img  :src="`data:image/png;base64,${imgB64}`" alt="分子结构图"  /><br>
                    </div>
                  </div>
                </el-card>


                



                


            </div>


      </div>
    </div>
  </div>
  
</template>

<script setup>
import { ref,watch } from 'vue';
import axios from 'axios';

const moleculeInput = ref('');  

const target_image_f = ref('')

const result_img = ref('')

const recentResults0 = ref([]); 
const recentResults1 = ref([]); 
const recentResults2 = ref([]); 

const scores_list0 = ref([])  // 分数
const scores_list1 = ref([])  // 分数
const scores_list2 = ref([])  // 分数



const first_edit = ref('')
const second_edit = ref('')


const inputs = ref([{ value: '' }])
const activeTab = ref("create")   // 表的选择


// 合成子预览的变量设置
const oneAtom = ref('')
const oneBond1 = ref('')
const oneBond2 = ref('')
const first_edit1 = ref('')
const first_edit2 = ref('')
const second_edit1 = ref('')
const second_edit2 = ref('')

const heCzi_img_danyuanzi = ref('')
const heCzi_img_danjian = ref('')
const heCzi_img = ref('')

const beam_num0 = ref()
const beam_num1 = ref()
const beam_num2 = ref()
const beam_num22 = ref()



//



async function TragetImage() {
  try {
    console.log('ok');
    const response = await axios.post('http://127.0.0.1:8000/user/1/TargetImage_canonicalize/',
        {
        // 发送的数据
        smile: String(moleculeInput.value)
        }
        
    );
    // console.log(molecule);  
    if (response.data && response.data.images) {
      target_image_f.value = response.data.images; // 保存最近结果
      console.log(response)
      // console.log(recentResults.value)
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

// 合成子预览
// 单原子
async function HeChengZi_oneAtmo() {
  console.log('合成子预-单原子')
  try{
      const response = await axios.post('http://127.0.0.1:8000/user/1/HeChengZi/',
        {
        // 发送的数据
        smile: String(moleculeInput.value),
        edit_1:String(oneAtom.value)+':0',
        edit_2:'0:0'
        }
      );

  if (response.data && response.data.img_heCzi) {
      heCzi_img_danyuanzi.value = response.data.img_heCzi; // 保存最近结果
      console.log(response)
    }
        
    

    
  }catch (error) {
    console.error('断键失败:', error);
  }
}
// 单键
async function HeChengZi_oneBond() {
  console.log('合成子预览-单键')
  try{
      const response = await axios.post('http://127.0.0.1:8000/user/1/HeChengZi/',
        {
        // 发送的数据
        smile: String(moleculeInput.value),
        edit_1:String(oneBond1.value)+':'+String(oneBond2.value),
        edit_2:'0:0'
        }
      );
      

  if (response.data && response.data.img_heCzi) {
      heCzi_img_danjian.value = response.data.img_heCzi; // 保存最近结果
      console.log(response)
    }
        
    

    
  }catch (error) {
    console.error('断键失败:', error);
  }
}
// 双
async function HeChengZi() {
  console.log('合成子预览')
  try{
      const response = await axios.post('http://127.0.0.1:8000/user/1/HeChengZi/',
        {
        // 发送的数据
        smile: String(moleculeInput.value),
        edit_1:String(first_edit1.value)+':'+String(first_edit2.value),
        edit_2:String(second_edit1.value)+':'+String(second_edit2.value),
        }
      );
      

  if (response.data && response.data.img_heCzi) {
      heCzi_img.value = response.data.img_heCzi; // 保存最近结果
      console.log(response)
    }
        
    

    
  }catch (error) {
    console.error('断键失败:', error);
  }
}


// 断键预测  任务提交
async function edit_context() {
  console.log('断键请求',activeTab.value);
  try{
      if (activeTab.value === 'danyuanzi'){
            const response = await axios.post('http://127.0.0.1:8000/user/1/edits_context002/',
            {
            // 发送的数据
            smile: String(moleculeInput.value),
            edit_1:String(oneAtom.value)+':0',
            beam_num:parseInt(beam_num0.value, 10)
            }
            );

            if (response.data && response.data.result_b64) {
                recentResults0.value = response.data.result_b64; // 保存最近结果
                scores_list0.value = response.data.scores_list;
                console.log(response)
              }
       }
      else if (activeTab.value === 'danjian'){
            const response = await axios.post('http://127.0.0.1:8000/user/1/edits_context002/',
                    {
                    // 发送的数据
                    smile: String(moleculeInput.value),
                    edit_1:String(oneBond1.value)+':'+String(oneBond2.value),
                    beam_num:parseInt(beam_num1.value, 10)
                    }
                  );
                  

            if (response.data && response.data.result_b64) {
                recentResults1.value = response.data.result_b64; // 保存最近结果
                scores_list1.value = response.data.scores_list;
                console.log(response)
              }
       }
      else if (activeTab.value === 'shuangjian'){
        console.log('双键预测执行')
              const response = await axios.post('http://127.0.0.1:8000/user/1/edits_context/',
                      {
                      // 发送的数据
                      smile: String(moleculeInput.value),
                      edit_1:String(first_edit1.value)+':'+String(first_edit2.value),
                      edit_2:String(second_edit1.value)+':'+String(second_edit2.value),
                      beam_num01:parseInt(beam_num2.value, 10),
                      beam_num02:parseInt(beam_num22.value, 10),
                      }
                    );
                    

            if (response.data && response.data.result_b64) {
                recentResults2.value = response.data.result_b64; // 保存最近结果
                scores_list2.value = response.data.scores_list;
                console.log(response)
              }
       }


    //   const response = await axios.post('http://127.0.0.1:8000/user/1/edits_context/',
    //     {
    //     // 发送的数据
    //     smile: String(moleculeInput.value),
    //     edit_1:String(first_edit.value),
    //     edit_2:String(second_edit.value)
    //     }
        
    // );
  
  // if (response.data && response.data.result_b64) {
  //     result_img.value = response.data.result_b64; // 保存最近结果
  //     console.log(response)
  //   }
    
  }catch (error) {
    console.error('断键失败:', error);
  }

}



const addInput = () => {
  inputs.value.push({ value: '' });
};

const removeInput = (index) => {
  inputs.value.splice(index, 1);
};


</script>

<style scoped>

.header {
  background-color: rgb(1,70,138);
  padding: 0px;
  text-align: center;
  font-size: 1.5em;  /* 文字大小 */
  font-weight: bold;  /* 加粗 */
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
  border: 2px solid #cccccc;  /* 设置 1 像素实线边框，颜色为浅灰色 */
  padding: 5px;
  flex: 1;
  display: flex;              /* 使用 Flexbox */
  flex-direction: column;     /* 垂直排列子元素 */
  
}

.submit-button {    /* 按键 */
  margin-top: 5px;
  margin-bottom: 10px;  
  /* line-height: 40px;  */
  /* padding: 5px 10px;  调整内边距，减少按钮的高度和宽度 */
  font-size: 0.6em;    /* 调整字体大小 */
  /* min-width: 100px;    可选: 设置最小宽度 */
  max-width: 100px;
  height: 50px;
  
}

.input_smiles{
    margin-top: 10px;
    margin-bottom: 10px;  

}
.recent_tarimage{
    max-width: 65%; /* 图像最大宽度为其父容器的 100% */
    height: auto;    /* 高度自适应，保持宽高比 */
    
    

}

.heCzi_image{
    max-width: 80%; /* 图像最大宽度为其父容器的 100% */
    height: auto;    /* 高度自适应，保持宽高比 */
    
    

}

/* 垂直布局，左对齐，顶部对齐 */
.card1-layout1{
  display: flex;
  flex-direction: column;
  align-items: flex-start; 
  justify-content: flex-start; 
  
  gap: 10px; 
}

.card1-layout2{
  display: flex; /* 启用 Flexbox 布局 */
  flex-direction: row; /* 水平排列子元素 */
  align-items: center; /* 垂直居中对齐子元素 */
  
  gap: 10px; /* 元素之间的间距 */
}

</style>