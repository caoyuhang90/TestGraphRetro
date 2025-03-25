<template>
  <el-container style="height: 100vh;">
    <el-header class="header">
      合成路线智能预测
    </el-header>
    <el-main style="display: flex; flex-direction: column; height: 100%;">

      <!-- SMILES 上传卡片 -->
      <el-card class="upload-card">
        <h3 class="title">上传信息</h3>
        <el-row :gutter="20" class="button-row">
          <el-button type="primary" @click="activeTab = 'smiles'">smiles反应式</el-button>
        </el-row>
        <el-form v-if="activeTab === 'smiles'" :model="form" label-width="100px" @submit.prevent="getPredictedYield">
          <el-form-item label="smiles">
            <el-input v-model="form.smiles" placeholder="输入化学反应的 SMILES 表达式"></el-input>
          </el-form-item>
          <el-form-item>
            <el-button type="primary" @click="getPredictedYield">提交</el-button>
          </el-form-item>
        </el-form>
      </el-card>

      <!-- 合成路线结果卡片 -->
      <el-card class="result-card" :style="{ flex: imagePaths.length === 0 ? '1' : 'unset' }">
        <h2 class="title" style="font-size: 32px !important;">合成路线</h2>
        <el-select v-if="imagePaths.length > 0" v-model="selectedImage" placeholder="选择结果图">
          <el-option v-for="(img, index) in imagePaths" :key="index" :label="'图像 ' + (index + 1)" :value="index"></el-option>
        </el-select>
        <div v-if="imagePaths.length > 0">
          <!-- 显示图像 -->
          <img :src="'http://localhost:5003/get_image/' + imagePaths[selectedImage]" alt="result image" style="max-width: 100%; height: auto;">
        </div>
        <p v-else class="result-title">请输入 产物SMILES 表达式并提交以查看预测结果。</p>
      </el-card>

    </el-main>
  </el-container>
</template>

<script>
import axios from 'axios';

export default {
  data() {
    return {
      form: {
        smiles: ''
      },
      imagePaths: [],  // 存储返回的图像文件名
      selectedImage: 0,
      activeTab: 'smiles'
    };
  },
  methods: {
    getPredictedYield() {
      // 清空之前的图像
      this.imagePaths = [];
      this.selectedImage = 0;

      axios.post('http://localhost:5003/generate_route', { smiles: this.form.smiles })
        .then(response => {
          if (response.data.image_paths) {
            this.imagePaths = response.data.image_paths;  // 存储图像文件名
          } else {
            this.$message.error('没有生成结果图像');
          }
        })
        .catch(error => {
          this.$message.error('请求失败: ' + error.message);
        });
    }
  }
};
</script>

<style scoped>
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

.result-card {
  margin-top: 20px;
  padding: 20px;
  text-align: center;
  border: 1px solid #dcdfe6;
  border-radius: 8px;
  flex: 1;
}


.result-card p {
  font-size: 24px;
  color: #5D9CEC;
  font-weight: bold;
}

.result-title {
  font-size: 32px !important;
  color:rgb(1, 70, 138) !important;
}
:deep(.el-button--primary) {
  background-color: rgb(1, 70, 138) !important;
  color: #ffffff !important;
  border-color: rgb(1, 70, 138) !important;
}
:deep(.el-form-item__label) {
  font-size: 20px;
}

/* 调整 el-select 和 el-option 的字体大小 */
:deep(.el-select .el-input__inner),
:deep(.el-option) {
  font-size: 20px;
}

/* 调整 el-option 的下拉菜单字体大小 */
:deep(.el-select-dropdown__item) {
  font-size: 20px;
}
</style>
