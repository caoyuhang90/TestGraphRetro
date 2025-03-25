<template>
  <el-container style="height: 100vh;">
    <el-header class="header">
      合成条件预测
    </el-header>
    <el-main style="display: flex; flex-direction: column; height: 100%;">
      <el-card class="upload-card">
        <h3>上传信息</h3>
        <el-row :gutter="20" class="button-row">
          <el-button type="primary" @click="activeTab = 'manual'">手动输入</el-button>
          <el-button type="info" disabled>图片/PDF</el-button>
          <el-button type="info" disabled>文本信息</el-button>
        </el-row>
        <el-form v-if="activeTab === 'manual'" :model="form" label-width="100px">
          <el-form-item label="smiles">
            <el-input v-model="form.smiles" placeholder="输入 SMILES 表达式，每行一个"></el-input>
          </el-form-item>
          <el-form-item>
            <el-button type="primary" @click="getPredictedConditions">提交</el-button>
          </el-form-item>
        </el-form>

        <el-upload
          class="upload-demo"
          action="http://localhost:5004/upload"
          :on-success="handleUploadSuccess"
          :before-upload="beforeUpload"
          :show-file-list="false"
          accept=".txt"
        >
          <template #trigger>
            <el-button type="primary">选择文件上传</el-button>
          </template>
        </el-upload>
      </el-card>

      <el-card class="result-card" :style="{ flex: results.length === 0 ? '1' : 'unset' }">
        <h3>预测的合成条件</h3>
        <div v-if="results.length === 0" class="no-results">
          <p>请输入数据以查看预测的合成条件。</p>
        </div>
        <el-table v-else :data="results" border style="width: 100%" class="styled-table">
          <el-table-column prop="index" label="top-k" width="60">
            <template #default="scope">
              {{ scope.$index + 1 }}
            </template>
          </el-table-column>
          <el-table-column prop="catalyst" label="Catalyst"></el-table-column>
          <el-table-column prop="solvent1" label="Solvent 1"></el-table-column>
          <el-table-column prop="solvent2" label="Solvent 2"></el-table-column>
          <el-table-column prop="reagent1" label="Reagent 1"></el-table-column>
          <el-table-column prop="reagent2" label="Reagent 2"></el-table-column>
          <el-table-column prop="score" label="Score">
            <template #default="scope">
              <span class="score-cell">{{ scope.row.score }}</span>
            </template>
          </el-table-column>
        </el-table>
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
      activeTab: 'manual',
      results: [],
    };
  },
  methods: {
    async getPredictedConditions() {
      try {
        const response = await axios.post('http://localhost:5004/predict', {
          smiles: this.form.smiles.split('\n'),
        });
        this.results = response.data;
      } catch (error) {
        console.error("Failed to fetch predicted conditions:", error);
        this.$message.error("无法获取预测条件，请重试");
      }
    },
    handleUploadSuccess(response) {
      this.form.smiles = response.smiles.join('\n');
      this.getPredictedConditions();
    },
    beforeUpload(file) {
      const isTxt = file.type === 'text/plain';
      if (!isTxt) {
        this.$message.error('只能上传 .txt 文件');
      }
      return isTxt;
    }
  },
};
</script>

<style scoped>
.header {
  background-color: rgb(1,70,138);
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
  border: 1px solid #dcdfe6;
  border-radius: 8px;
  box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
  background-color: white;
  flex: 1; /* 使结果卡片填满剩余空间 */
}

.no-results {
  text-align: center;
  color: #999;
}

.styled-table th {
  background-color: #4CAF50;
  color: white;
  text-align: left;
}

.styled-table tr:nth-child(even) {
  background-color: #f9f9f9;
}

.styled-table tr:hover {
  background-color: #f1f1f1;
}

.styled-table th, .styled-table td {
  padding: 12px;
  border: 1px solid #ddd;
  text-align: left;
}

.score-cell {
  font-weight: bold;
  color: #d9534f;
}
</style>
