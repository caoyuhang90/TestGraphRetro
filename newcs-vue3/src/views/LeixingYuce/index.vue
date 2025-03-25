<template>
  <el-container style="height: 100vh;">  <!-- 全部子元素会垂直上下排列 -->
    <!-- 顶部导航栏 -->
    <el-header class="header">
      化学反应类型预测
    </el-header>

    <!-- 主体内容 -->
    <!-- 输入方式模块 -->
    <div class="container" v-if="activeModule === 'input'">
      <el-card class="module" style="width: 100%; margin-right: 30px; height: 900px; position: relative;">
        <h2 class="method-title">选择输入方式</h2>
        <el-radio-group v-model="inputMethod">  <!-- radio 单选框  选择模式 -->
          <el-radio class="radio" :value="'string'">输入SMILES字符串</el-radio>  <!-- :value会传给v-model的值 -->
          <el-radio class="radio" :value="'file'">上传SMILES文件</el-radio>
        </el-radio-group>

        <!-- Input Methods -->
        <el-form label-position="top" style="margin-top: 20px;">  <!-- 创建表单，用于收集数据 -->
          <el-form-item label="选择输入方式:">  
            <el-select v-model="inputMethod" placeholder="选择输入方式">  <!-- 下拉框 -->
              <el-option label="输入SMILES字符串" value="string"></el-option>
              <el-option label="输入反应物、反应条件和产物" value="reaction"></el-option>
            </el-select>
          </el-form-item>

          <!-- SMILES String Input -->
          <el-form v-if="inputMethod === 'string'" label-position="top">  <!-- 表单中，负责包裹单个表单控件（如输入框、选择框等）及其标签 -->
            <el-form-item label="输入SMILES字符串:">
              <el-input v-model="smilesString" placeholder="请输入SMILES字符串"></el-input>
            </el-form-item>
          </el-form>

          <!-- Reaction Inputs -->
          <el-form v-if="inputMethod === 'reaction'" label-position="top">
            <el-form-item label="输入反应物:">
              <el-input v-model="reactants" placeholder="请输入反应物"></el-input>
            </el-form-item>
            <el-form-item label="输入反应条件:">
              <el-input v-model="reactionConditions" placeholder="请输入反应条件"></el-input>
            </el-form-item>
            <el-form-item label="输入产物:">
              <el-input v-model="products" placeholder="请输入产物"></el-input>
            </el-form-item>
          </el-form>

          <!-- Model Upload -->
          <el-form-item label="上传预测模型（选填，默认使用默认模型）:">
            <el-upload class="upload-demo" 
            drag  
            :before-upload="() => false"
            @change="handleFileChange($event, 'modelFile')">
              <i class="el-icon-upload"></i>
              <div class="el-upload__text">点击或拖拽文件上传</div>  
            </el-upload>
          </el-form-item>

          <!-- Submit Button -->
          <el-form-item>
            <el-button class="submit-button" type="primary" @click="submitPrediction">提交信息</el-button>
          </el-form-item>
        </el-form>

        <!-- 切换按钮，位于左下角 如果activeModule='input'就返回'训练模型'，否则返回'输入方式' -->
        <el-button class="toggle-button" type="primary" @click="toggleModule" style="position: absolute; bottom: 10px; left: 10px;">
          切换到 {{ activeModule === 'input' ? '训练模型' : '输入方式' }}
        </el-button>
      </el-card>

      <!-- Prediction Result Card -->
      <el-card class="result-card" style="width: 100%; margin-right: 30px; height: 940px;">
        <h3>预测结果</h3>
        <p v-if="predictionResult"><strong>Class Name:</strong> {{ predictionResult.class_name }}</p>
      </el-card>
    </div>

    <!-- 添加新反应类型数据训练模型模块 -->
    <el-card v-if="activeModule === 'training'" class="module" style="width: 100%; margin-top: 0px; height: 900px; position: relative;">
      <h2>添加新反应类型数据训练模型</h2>
      <el-form label-position="top" ref="form" :model="form">
        <el-form-item label="上传新反应类型数据:">
          <el-upload class="upload-demo" 
          drag 
          :before-upload="() => false"
          @change="handleFileChange($event, 'newRxnFile')">
            <i class="el-icon-upload"></i>
            <div class="el-upload__text">点击或拖拽文件上传</div>
          </el-upload>
        </el-form-item>

        <el-form-item label="反应类别 ID:">
          <el-input v-model="form.classId" placeholder="输入反应类别 ID"></el-input>
        </el-form-item>
        <el-form-item label="反应类别编号:">
          <el-input v-model="form.rxn_class" placeholder="输入反应类别编号"></el-input>
        </el-form-item>
        <el-form-item label="反应类别名称:">
          <el-input v-model="form.className" placeholder="输入反应类别名称"></el-input>
        </el-form-item>

        <el-form-item>
          <el-button class="submit-button" type="primary" @click="submitTraining">开始训练</el-button>
        </el-form-item>
      </el-form>

      <!-- 切换按钮，位于左下角 -->
      <el-button class="toggle-button" type="primary" @click="toggleModule" style="position: absolute; bottom: 10px; left: 10px;">
        切换到 {{ activeModule === 'input' ? '训练模型' : '输入方式' }}
      </el-button>
    </el-card>
  </el-container>
</template>



<script >
import { Compass } from '@element-plus/icons-vue';
import { ElMessage } from 'element-plus';
import axios from 'axios';

export default {
  name: 'SynthesisPrediction',
  components: {
    Compass
  },
  data() {
    return {
      activeModule: 'input', // 默认为选择输入方式模块
      inputMethod: 'string', // 默认为输入字符串的方式
      smilesString: '',
      reactants: '',
      reactionConditions: '',
      products: '',
      modelFile: null,
      smilesFile: null,
      modelFile2: null,
      newRxnFile: null,
      taskName: null,
      form: {  // 初始化 form 对象
        taskName: '',
        classId: '',  // 初始化 classId
        className: '',  // 初始化 className
        rxn_class: ''
      },
      predictionResult: null,  // 存储预测结果
      fileData: '' // 存储文件数据
    };
  },

  methods: {
    toggleModule() {
      // 切换模块显示
      this.activeModule = this.activeModule === 'input' ? 'training' : 'input';
    },
    handleFileChange(file, fieldName) {
      if (file) {
        this[fieldName] = file.raw || file;  // 使用 file.raw 或 file
      } else {
        ElMessage.error('请上传文件');
      }
    },

    async submitPrediction() {
      if (this.inputMethod === 'reaction') {
        if (!this.reactants || !this.reactionConditions || !this.products) {
          ElMessage.error('请提供完整的反应物、反应条件和产物信息');
          return;
        }
        this.smilesString = `${this.reactants}.${this.reactionConditions}>>${this.products}`;
      }
      if (!this.smilesString && !this.modelFile) {
        ElMessage.error('请提供SMILES字符串或上传模型文件');
        return;
      }
      const formData = new FormData();
      formData.append('smiles_string', this.smilesString);
      if (this.modelFile) formData.append('model_file', this.modelFile);

      try {
        const response = await axios.post('http://localhost:5001/predict', formData);
        this.predictionResult = response.data;  // 将结果存储在 predictionResult 中
        ElMessage.success('预测成功');
      } catch (error) {
        ElMessage.error('预测失败');
        console.error('Error:', error);
      }
    },

    async submitFilePrediction() {
      if (!this.smilesFile) {
        ElMessage.error('请上传SMILES文件');
        return;
      }

      const formData = new FormData();
      formData.append('smiles_file', this.smilesFile);
      if (this.modelFile2) {
        formData.append('model_file', this.modelFile2);
      }

      try {
        const response = await axios.post('http://localhost:5001/file_predict', formData, { responseType: 'blob' });
        const blob = new Blob([response.data], { type: 'text/csv' });
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'predicted_results.csv';
        document.body.appendChild(a);
        a.click();
        window.URL.revokeObjectURL(url);
        ElMessage.success('文件预测成功');
      } catch (error) {
        ElMessage.error('文件预测失败');
        console.error('Error:', error);
      }
    },

    async submitTraining() {
      if (!this.newRxnFile || !this.form.classId || !this.form.className) {
        ElMessage.error('请提供任务名称、反应类别ID、反应类别名称和上传数据');
        return;
      }

      const formData = new FormData();
      formData.append('new_rxn_file', this.newRxnFile);
      formData.append('class_id', this.form.classId);
      formData.append('class_name', this.form.className);
      formData.append('task_name2', this.form.taskName);
      formData.append('rxn_class', this.form.rxn_class);

      try {
        const response = await axios.post('http://localhost:5001/train', formData, {
          responseType: 'blob'
        });

        const url = window.URL.createObjectURL(new Blob([response.data]));
        const link = document.createElement('a');
        link.href = url;
        link.setAttribute('download', 'model_files.zip');
        document.body.appendChild(link);
        link.click();
        window.URL.revokeObjectURL(url);
        document.body.removeChild(link);
        ElMessage.success('训练任务成功提交，ZIP 文件已下载');
      } catch (error) {
        ElMessage.error('训练任务提交失败');
        console.error('Error:', error.response);
      }
    }
  },
};
</script>

<style scoped>
/*  */
.header {
  background-color: rgb(1, 70, 138);   /* 背景颜色 */
  color: white;  /* 字体颜色 */
  text-align: center;  /* text-align -- 文字排版    居中   */
  padding: 0px 0;   /* 与其他边框的内边距 */
}

.el-tabs {
  background-color: #ffffff;
  margin-bottom: 20px;
}

.module {
  padding: 20px;
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
  border-radius: 8px;  /* 边框圆角 */
}

.el-upload__text {
  font-size: 20px;
  /* 点击或拖拽文件上传 */
  text-align: center;
}

.el-button.toggle-button {
  width: 150px;
  height: 50px;
  margin-bottom: 0;
  margin-left: 20px;
  background-color: #3d5b88;
  /* 淡绿色背景 */
  text-align: left;
  /* 左对齐 */
  border: none;
  /* 无边框 */
  font-size: 18px;
}

.submit-button {
  color: rgb(1, 70, 138);

}

.reaction-tooltip .reaction-list {
  background-color: #ffffff;
  /* 白色背景 */
  color: #000000;
  /* 黑色字体 */
  max-height: 200px;
  /* 控制最大高度，约10行 */
  overflow-y: auto;
  /* 允许垂直滚动 */
  width: 250px;
  /* 可根据需求调整宽度 */
  border: 1px solid #dcdcdc;
  /* 边框颜色，可选 */
  border-radius: 4px;
  /* 圆角效果，可选 */
  box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.1);
  /* 轻微阴影，可选 */
}

.reaction-tooltip .reaction-list ul {
  list-style: none;
  /* 移除默认的列表样式 */
  padding: 0;
  margin: 0;
}

.reaction-tooltip .reaction-list li {
  padding: 5px 10px;
  /* 控制每项的间距 */
  font-size: 14px;
  /* 控制字体大小 */
  line-height: 1.5;
}

.result-card {
  font-size: 32px;
}

/* 调整选择输入方式标题的【字体大小】 */
.method-title {
  font-size: 32px;  
  /* 根据需要调整 */
}

/* 调整 el-radio 的标签文字字体大小 */
:deep(.el-radio__label) {
  font-size: 20px;
  /* 输入SMILES字符串 */
}

/* 调整 el-form-item 标签的字体大小 */
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

.container {
  display: flex;  /* 启用 Flexbox 布局 */
  flex-direction: row;  
  flex-wrap: nowrap;
  width: 100%;
}
</style>

