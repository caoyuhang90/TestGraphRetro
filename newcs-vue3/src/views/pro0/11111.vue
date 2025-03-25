<template>
  <div>
    <el-header class="header">
      合成路线智能预测
    </el-header>

    <el-tabs v-model="activeTab" class="tab-container">
      <el-tab-pane label="创建任务" name="create">
        <div class="container">
          <div class="column">
            <!-- 上传信息模块 -->
            <div class="section upload-section">
              <h3>上传信息</h3>
              <el-input v-model="moleculeInput" placeholder="输入SMILES字符串" clearable style="width: 52%;"></el-input>
              <el-input v-model="number_of_routes" placeholder="输入期望的路线数量" clearable style="width: 52%;"></el-input>
              <el-button type="primary" class="submit-button" @click="submitTask">开始预测</el-button>
            </div>

            <!-- 预测结果模块 -->
            <div class="section prediction-results" v-if="predictionResults">
              <h3>预测结果</h3>
              <hr>
              <div v-for="(imgURL, index) in list" :key="index" class="reaction-item">
                <div class="reaction-info">
                  <h4>index:{{ index }}</h4>
                  <img :src="imgURL" alt="分子结构图" class="small-image" @click="toggleExpand(index)"/><br>
                </div>

                <el-button @click="toggleExpand(index)">
                  {{ expandedIndex.value === index ? '折叠' : '展开' }}反应底物大图
                </el-button>
                
                <div v-if="expandedIndex.value === index" class="overlay">
                  <img :src="imgURL" alt="大图" class="large-image" />
                  <el-button @click="toggleExpand(index)">关闭</el-button>
                </div>
              </div>
              <hr>
            </div>
          </div>
        </div>
      </el-tab-pane>

      <el-tab-pane label="最近结果" name="recent">
        <!-- 最近结果的内容 -->
      </el-tab-pane>
    </el-tabs>
  </div>
</template>

<script setup>
import { ref } from 'vue';
import axios from 'axios';

const moleculeInput = ref('');
const number_of_routes = ref(1);
const totalImages_f = ref(0);
const top_score_f = ref(0);
const activeTab = ref("create");
const predictionResults = ref(null);
const expandedIndex = ref(null);
const result = ref();
const list = ref([]);

async function submitTask() {
    const response = await axios.post('http://127.0.0.1:3000/api/generate-molecule', {
        molecule: { molecule: String(moleculeInput.value) },
        numb_routes: number_of_routes.value,
    });

    totalImages_f.value = response.data.total_images;
    top_score_f.value = response.data.top_score;

    if (response.data.images && response.data.images.length > 0) {
        list.value = response.data.images.map(img => `data:image/png;base64,${img}`);
    }
    predictionResults.value = true;
}

function toggleExpand(index) {
    expandedIndex.value = (expandedIndex.value === index) ? null : index;
}
</script>

<style>
.reaction-item {
    position: relative; /* 使大图能够覆盖 */
}

.overlay {
    position: absolute;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    background-color: rgba(255, 255, 255, 0.9); /* 半透明背景 */
    display: flex;
    justify-content: center;
    align-items: center;
}

.large-image {
    max-width: 100%;
    max-height: 100%;
}
</style>