<template>
  <el-container style="height: 100vh;">
    <el-header class="header">
      合成反应产率预测
    </el-header>
    <el-main style="display: flex; flex-direction: column; height: 100%;">
      <el-card class="upload-card">
        <h3>上传信息</h3>
        <el-row :gutter="20" class="button-row">
          <el-button type="primary" @click="activeTab = 'smiles'">smiles反应式</el-button>
          <el-button type="info" disabled>图片/PDF</el-button>
          <el-button type="info" disabled>文本信息</el-button>
        </el-row>
        <el-form v-if="activeTab === 'smiles'" :model="form" label-width="100px" @submit.prevent="getPredictedYield">
          <el-form-item label="smiles表达式">
            <el-input v-model="form.smiles" placeholder="输入化学反应的 SMILES 表达式"></el-input>
          </el-form-item>
          <el-form-item>
            <el-button type="primary" @click="getPredictedYield">提交</el-button>
          </el-form-item>
        </el-form>
      </el-card>

      <el-card class="result-card" :style="{ flex: predictedYield === null ? '1' : 'unset' }">
        <h3>预测的反应产率</h3>
        <p v-if="predictedYield !== null">{{ predictedYield }}%</p>
        <p v-else>请输入 SMILES 表达式并提交以查看预测结果。</p>
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
        smiles: '',
      },
      activeTab: 'smiles',
      predictedYield: null,
    };
  },
  methods: {
    async getPredictedYield() {
      try {
        const response = await axios.post('http://localhost:5005/predict', {
          smiles: this.form.smiles,
        });
        this.predictedYield = response.data.predicted_yield;
      } catch (error) {
        console.error("Failed to fetch predicted yield:", error);
        this.$message.error("无法获取预测产率，请重试");
      }
    },
  },
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

.el-button {
  background-color: #C1E7E3;
  color: #333;
}

.result-card {
  margin-top: 20px;
  padding: 20px;
  text-align: center;
  border: 1px solid #dcdfe6;
  border-radius: 8px;
  flex: 1; /* 使结果卡片填满剩余空间 */
}

.result-card h3 {
  font-size: 24px;
  margin-bottom: 10px;
}

.result-card p {
  font-size: 24px;
  color: #5D9CEC;
  font-weight: bold;
}
</style>
